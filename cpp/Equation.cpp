
/*
------------------------------------------------------------------------------------------------------
	Файл:		Equation.cpp
	Версия:		1.06
	DLM:		10.04.2004
------------------------------------------------------------------------------------------------------
*/

#include <fstream.h>
#include <math.h>
#include <assert.h>

#include "Equation.h"
#include "Engine.h"

void SolveODU(CInterval AB, double U0, double (*fn)(double, double), char *fname)
{
	double
		*U = new double[AB.N()],
		h = AB.H();
	assert(U!=0);
	U[0] = U0;

	double k1,k2,k3,k4,x;
	for(int i=0; i<AB.N()-1; i++)
	{
		x = AB.X(i);
		k1 = fn(x, U[i]);
		k2 = fn(x + h/2., U[i] + h*k1/2.);
		k3 = fn(x + h/2., U[i] + h*k2/2.);
		k4 = fn(x + h, U[i] + h*k3);
		U[i+1] = U[i] + h/6.*( k1+2.*( k2+k3 )+k4 );
	}

	ofstream f(fname);
	for(i=0; i<AB.N(); i++) f << AB.X(i) <<" "<< U[i] <<"\n";
	f.close();

	delete []U;
}

CEquation::CEquation(double x1, double x2, double t2, int l, int r, double lH, double rH)
{
	assert(l>=1 && l<=3 && r>=1 && r<=3 && lH>0 && rH>0);

	D.x = new CInterval(x1, x2, 1);
	D.t = new CInterval(0, t2, 1);
	assert(D.x!=0 && D.t!=0);
	
	flag = false; lbt = l; rbt = r; lh = lH; rh = rH;
}

CEquation::~CEquation()
{
	if(flag) delete Arr;
	delete D.x;
	delete D.t;
}

void CEquation::Solve(double h, double tau)
{
	//установка шагов по пространственной и временной координатам
	assert(h>0 && tau>0);
	D.x->ReBorn(D.x->X1(), D.x->X2(), h);
	D.x->ReBorn(D.x->X1(), D.x->X2(), tau);
	
	//выделение места под решение уравнения
	if(flag) delete Arr;
	Arr = new CMatrix(D.x->N()+1, D.x->N()+1);
	assert(Arr);

	flag = true;

	//выделение места под коэффициенты для метода прогонки
	int N = D.x->N();
	double
		*A = new double[N],	*B = new double[N],
		*C = new double[N],	*F = new double[N],
		*U = new double[N+1];
	assert(A!=0 && B!=0 && C!=0 && F!=0 && U!=0);
	
	//задание начального значения
	for(int i=0; i<=N; i++) Arr->get(0,i) = gU0(i);
	
	//
	for(int j=0; j<D.x->N(); j++)
	{
		for(i=1; i<N; i++)
		{
			double
				_u = Arr->get(j,i-1),
				u = Arr->get(j,i),
				u_ = Arr->get(j,i+1),

				Alpha =	( gK(i,j,u) + gK(i+1,j,u_) + gV(i,j,u)*h )/2.,
				Beta =	( gK(i,j,u) + gK(i-1,j,_u) - gV(i,j,u)*h )/2.,
				Gamma = 2.*pow(h,2)/tau*gM(i,j,u);

			A[i] = Beta;
			B[i] = Alpha;
			C[i] = Alpha + Beta + Gamma;
			F[i] = u_*Alpha + _u*Beta - u*(Alpha + Beta - Gamma) + 2.*pow(h,2)*gF(i,j,u);
		}

		double
			Mu[3] = {0},
			Nu[3] = {0};
		int nj = j+1;
		switch(lbt)
		{
			case 1:	Nu[1] = gLU(nj); break;
			case 2:	Mu[1] = 1.; Nu[1] = -h*gLdU_dx(nj); break;
			case 3: Mu[1] = 1./(1. + h*lh); Nu[1] = h*lh*gLTeta(nj)/(1. + h*lh); break;
		}
		switch(lbt)
		{
			case 1:	Nu[2] = gRU(nj); break;
			case 2:	Mu[2] = 1.; Nu[2] = h*gRdU_dx(nj); break;
			case 3: Mu[2] = 1./(1. - h*rh); Nu[2] = -h*rh*gRTeta(nj)/(1. - h*rh); break;
		}

		Progonka(N, A, B, C, F, Mu[1], Nu[1], Mu[2], Nu[2], U);
		for(i=0; i<=N; i++) Arr->get(nj,i) = U[i];
	}
	
	delete []A; delete []B; delete []C; delete []F;	delete []U;
}

void CEquation::sUt(char *fname, double t[], int size)
{
	assert(flag);
	for(int i=0; i<size; i++) assert(t[i]>=D.x->X1() && t[i]<=D.x->X2());

	ofstream f(fname);
	for(i=0; i<=D.x->N(); i++)
	{
		f<< D.x->X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr->get(D.x->i(t[j]),i);
		f<<"\n";
	}
	f.close();
}

void CEquation::sUx(char *fname, double x[], int size)
{
	assert(flag);
	for(int i=0; i<size; i++) assert(x[i]>=D.x->X1() && x[i]<=D.x->X2());

	ofstream f(fname);
	for(i=0; i<=D.x->N(); i++)
	{
		f<< D.x->X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr->get(i,D.x->i(x[j]));
		f<<"\n";
	}
	f.close();
}

CMatrix* CEquation::gUt(int it)
{
	assert( it>=0 && it<Arr->GetN()	);
	CMatrix *M = new CMatrix(2,Arr->GetN());
	assert(M);

	for(int i=0; i<Arr->GetN(); i++)
	{
		M->get(0,i) = D.x->X(i);
		M->get(1,i) = Arr->get(it,i);
	}
	return M;
}

CMatrix* CEquation::gUx(int ix)
{
	assert( ix>=0 && ix<Arr->GetN() );
	CMatrix *M = new CMatrix(2,Arr->GetM());
	assert(M);

	for(int i=0; i<Arr->GetM(); i++)
	{
		M->get(0,i) = D.x->X(i);
		M->get(1,i) = Arr->get(i,ix);
	}
	return M;
}