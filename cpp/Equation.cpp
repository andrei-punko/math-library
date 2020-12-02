
/*
------------------------------------------------------------------------------------------------------
	Файл:		Equation.cpp
	Версия:		1.10
	DLM:		16.03.2005
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
		x = AB.X(0),
		U = U0,
		h = AB.H(),
		k1,k2,k3,k4;

	ofstream f(fname); f<< x <<" "<< U <<"\n";

	for(int i=0; i<AB.N()-1; i++)
	{
		x = AB.X(i);
		k1 = fn(x, U);
		k2 = fn(x + h/2., U + h*k1/2.);
		k3 = fn(x + h/2., U + h*k2/2.);
		k4 = fn(x + h, U + h*k3);

		f<< x <<" "<< (U += h/6.*( k1+2*(k2+k3)+k4 )) <<"\n";
	}
	f.close();
}

void SolveODU(CInterval AB, double *U0, double (*fn)(double, double *, int), int N, char *fname)
{
	double
		*tU = new double[N],
		*dU = new double[N],
		*U = new double[N],
		h = AB.H(),x;
	
	assert(tU!=0 && dU!=0 && U!=0);

	ofstream f(fname); f<< AB.X(0);
	for(int i=0; i<N; i++) f<<" "<< (U[i] = U0[i]);

	f<<"\n";

	for(i=0; i<AB.N()-1; i++)
	{
		f<< (x = AB.X(i));
		for(int k=0; k<N; k++) dU[k] = h*fn(x,U,k);
		for(int j=0; j<N; j++) f<<" "<< (U[j] += dU[j]);

		f<<"\n";
	}
	f.close();

	delete []tU; delete []dU; delete []U;
}

CEqn::CEqn(double x1, double x2, double t2, int l, int r, double lH, double rH)
{
	assert(l>=1 && l<=3 && r>=1 && r<=3 && lH>0 && rH>0);

	D.x = new CInterval(x1, x2, 1);
	D.t = new CInterval(0, t2, 1);
	assert(D.x!=0 && D.t!=0);
	
	lbt = l; rbt = r; lh = lH; rh = rH;
}

CEqn::~CEqn()
{
	delete D.x;	delete D.t;
}


void CEqn::sUt(char *fname, double t[], int size)
{
	for(int i=0; i<size; i++) assert(t[i]>=D.t->X1() && t[i]<=D.t->X2());

	ofstream f(fname);
	for(i=0; i<=D.x->N(); i++)
	{
		f<< D.x->X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr(D.t->i(t[j]),i);
		f<<"\n";
	}
	f.close();
}

void CEqn::sUx(char *fname, double x[], int size)
{
	for(int i=0; i<size; i++) assert(x[i]>=D.x->X1() && x[i]<=D.x->X2());

	ofstream f(fname);
	for(i=0; i<=D.t->N(); i++)
	{
		f<< D.t->X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr(i,D.x->i(x[j]));
		f<<"\n";
	}
	f.close();
}

void CEqn::gUt(int it, CMatrix &M)
{
	int N = Arr.GetN();
	assert(it>=0 && it<N);
	M.SetSize(2,N);

	for(int i=0; i<N; i++)
	{
		M(0,i) = D.x->X(i);
		M(1,i) = Arr(it,i);
	}
}

void CEqn::gUx(int ix, CMatrix &M)
{
	int N = Arr.GetM();	
	assert(ix>=0 && ix<N);
	M.SetSize(2,N);

	for(int i=0; i<N; i++)
	{
		M(0,i) = D.x->X(i);
		M(1,i) = Arr(i,ix);
	}
}

void CEqn::Solve(double h, double tau)
{
	//установка шагов по пространственной и временной координатам
	assert(h>0 && tau>0);
	D.x->ReBorn(D.x->X1(), D.x->X2(), h);
	D.t->ReBorn(D.t->X1(), D.t->X2(), tau);
	
	//выделение места под решение уравнения
	Arr.SetSize(D.t->N()+1, D.x->N()+1);

	//задание начального значения
	for(int i=0; i<=D.x->N(); i++) Arr(0,i) = gU0(i);
}

void CParEqn::Solve(double h, double tau)
{
	CEqn::Solve(h,tau);

	//выделение места под коэффициенты для метода прогонки
	int N = D.x->N();
	double
		*A = new double[N], *B = new double[N],
		*C = new double[N], *F = new double[N],
		*U = new double[N+1];
	assert(A!=0 && B!=0 && C!=0 && F!=0 && U!=0);
	
	//
	for(int j=0; j<D.t->N(); j++)
	{
		for(int i=1; i<N; i++)
		{
			double
				_u = Arr(j,i-1),
				u = Arr(j,i),
				u_ = Arr(j,i+1),

				Alpha =	( gK(i,j,u) + gK(i+1,j,u_) + gV(i,j,u)*h )/2.,
				Beta =	( gK(i,j,u) + gK(i-1,j,_u) - gV(i,j,u)*h )/2.,
				Gamma = 2.*pow(h,2)/tau*gM(i,j,u);

			A[i] = Beta;
			B[i] = Alpha;
			C[i] = Alpha + Beta + Gamma;
			F[i] = u_*Alpha + _u*Beta - u*(Alpha + Beta - Gamma) + 2.*pow(h,2)*gF(i,j,u);
		}
		
		int nj = j+1;
		double
			Mu[3] = {0},
			Nu[3] = {0};

		switch(lbt)
		{
			case 1:	Nu[1] = gLU(nj); break;
			case 2:	Mu[1] = 1.; Nu[1] = -h*gLdU_dx(nj); break;
			case 3: Mu[1] = 1./(1. + h*lh); Nu[1] = h*lh*gLTeta(nj)/(1. + h*lh); break;
		}
		switch(rbt)
		{
			case 1:	Nu[2] = gRU(nj); break;
			case 2:	Mu[2] = 1.; Nu[2] = h*gRdU_dx(nj); break;
			case 3: Mu[2] = 1./(1. - h*rh); Nu[2] = -h*rh*gRTeta(nj)/(1. - h*rh); break;
		}

		Progonka(N, A, B, C, F, Mu[1], Nu[1], Mu[2], Nu[2], U);
		for(i=0; i<=N; i++) Arr(nj,i) = U[i];
	}
	
	delete []A; delete []B; delete []C; delete []F;	delete []U;
}

//	CHypEquation
//
void CHypEqn::Solve(double h, double tau)
{
	CEqn::Solve(h,tau);
	
	//выделение места под коэффициенты для метода прогонки
	int N = D.x->N();
	double
		*A = new double[N],	*B = new double[N],
		*C = new double[N],	*F = new double[N],
		*U = new double[N+1];
	assert(A!=0 && B!=0 && C!=0 && F!=0 && U!=0);

	//задание граничных значений на первом слое
	Arr(1,0) = gLU(1);
	Arr(1,N) = gRU(1);

	//вычисление значения функции на первом слое для запуска разностной схемы
	for(int i=1; i<N; i++)
	{
		double
			_u = Arr(0,i-1),
			u = Arr(0,i),
			u_ = Arr(0,i+1);

		Arr(1,i) = u + tau*(gdU_dt(i) +
							tau/(2.*gM(i,0,u))*( 
							gK(i,0,u)/pow(h,2.)*(_u - 2.*u + u_) + gV(i,0,u)/(2.*h)*(u_ - _u) + gF(i,0,u) ));
	}

	//реализация разностной схемы
	for(int j=0; j<=D.t->N()-2; j++)
	{
		for(i=1; i<N; i++)
		{
			double
				_u = Arr(j,i-1),
				u = Arr(j,i),
				u_ = Arr(j,i+1),

				Alpha =	gK(i,j,u) - gV(i,j,u)*h/2.,
				Beta =	gK(i,j,u) + gV(i,j,u)*h/2.,
				Gamma = gL(i,j,u)*pow(h,2)/tau,
				Delta = 2*gM(i,j,u)*pow(h/tau,2);

			A[i] = Alpha;
			B[i] = Beta;
			C[i] = Alpha + Beta - Gamma + Delta;
			F[i] = _u*Alpha + u_*Beta - u*(Alpha + Beta + Gamma + Delta) +
					2*( Arr(j+1,i)*Delta + gF(i,j,u)*pow(h,2) );
		}

		int nj = j+2;
		double
			Mu[3] = {0},
			Nu[3] = {0};

		switch(lbt)
		{
			case 1:	Nu[1] = gLU(nj); break;
			case 2:	Mu[1] = 1; Nu[1] = -h*gLdU_dx(nj); break;
			case 3: Mu[1] = 1/(1 + h*lh); Nu[1] = h*lh*gLTeta(nj)/(1 + h*lh); break;
		}
		switch(rbt)
		{
			case 1:	Nu[2] = gRU(nj); break;
			case 2:	Mu[2] = 1; Nu[2] = h*gRdU_dx(nj); break;
			case 3: Mu[2] = 1/(1 - h*rh); Nu[2] = -h*rh*gRTeta(nj)/(1 - h*rh); break;
		}

		Progonka(N, A, B, C, F, Mu[1], Nu[1], Mu[2], Nu[2], U);
		for(i=0; i<=N; i++) Arr(nj,i) = U[i];
	}

	delete []A; delete []B; delete []C; delete []F;	delete []U;
}