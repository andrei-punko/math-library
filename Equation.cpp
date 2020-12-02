
/*
----------------------------------------------------------------------------------------
	Файл:		Equation.cpp
	Версия:		1.02
	DLM:		10.11.2003
----------------------------------------------------------------------------------------
*/

#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "Equation.h"
#include "Engine.h"
#include "Matrix.h"

void SolveODU(CInterval AB, const double U0, double (*fn)(double, double), const char *fname)
{
	double
		*U = new double[AB.N()],
		h = AB.H();
	U[0] = U0;

	double k1,k2,k3,k4,x;
	for(int i=0; i<AB.N()-1; i++)
	{
		x = AB.X(i);
		k1 = fn(x, U[i]);
		k2 = fn(x + h/2.0, U[i] + h*k1/2.0);
		k3 = fn(x + h/2.0, U[i] + h*k2/2.0);
		k4 = fn(x + h, U[i] + h*k3);
		U[i+1] = U[i] + h/6.0*( k1+2.0*( k2+k3 )+k4 );
	}

	ofstream f(fname);
	for(i=0; i<AB.N(); i++) f << AB.X(i) <<" "<< U[i] <<"\n";
	f.close();

	delete []U;
}

CEquation::CEquation(const double x1, const double x2, const double t2,
					 const int l, const int r, const double lH, const double rH)
{
	if(l<1 || l>3 || r<1 || r>3 || !(lH>0 && rH>0))
	{
		cout << "CEquation::CEquation(): error\n";
		exit(EXIT_FAILURE);
	}

	x12 = new CInterval(x1, x2, 1);
	t12 = new CInterval(0, t2, 1);
	
	flag = false; lbt = l; rbt = r; lh = lH; rh = rH;
}

CEquation::~CEquation()
{
	if(flag) delete Arr;
	delete x12;
	delete t12;
}

void CEquation::Solve(const double h, const double tau)
{
	x12->ReBorn(x12->X1(), x12->X2(), h);
	t12->ReBorn(t12->X1(), t12->X2(), tau);
	
	if(flag) delete Arr;
	Arr = new CMatrix<double>(t12->N()+1, x12->N()+1); flag = true;

	int N = x12->N();
	double
		*A = new double[N],	*B = new double[N],
		*C = new double[N],	*F = new double[N],
		*U = new double[N+1];
	
	for(int i=1; i<N; i++) Arr->get(0,i) = GetU0(x12->X(i));
	Arr->get(0,0) = GetLeftU(0);
	Arr->get(0,N) = GetRightU(0);
	
	for(int j=0; j< t12->N(); j++)
	{
		for(i=1; i<N; i++)
		{
			double
				_x = x12->X(i),
				_t = t12->X(j),
				_u = Arr->get(j,i),

		Alpha =	( GetK(_x,_t,_u) + GetK(x12->X(i+1),_t,Arr->get(j,i+1)) + GetV(_x,_t,_u)*h )/2,
		Beta =	( GetK(x12->X(i-1),_t,Arr->get(j,i-1)) + GetK(_x,_t,_u) - GetV(_x,_t,_u)*h )/2,
		Gamma = 2*pow(h,2)/tau*GetM(_x,_t,_u);

			A[i] = Beta;
			B[i] = Alpha;
			C[i] = Alpha + Beta + Gamma;
			F[i] = Arr->get(j,i+1)*Alpha + Arr->get(j,i-1)*Beta -
				_u*(Alpha + Beta - Gamma) + 2*pow(h,2)*GetF(_x,_t,_u);
		}

		double
			Mu[3] = {0},
			Nu[3] = {0},
			t = t12->X(j+1);
		switch(lbt)
		{
			case 1:	Nu[1] = GetLeftU(t); break;
			case 2:	Mu[1] = 1; Nu[1] = -h*GetLeftU_x(t); break;
			case 3: Mu[1] = 1/(1 + h*lh); Nu[1] = h*lh*GetLeftTeta(t)/(1 + h*lh); break;
		}
		switch(lbt)
		{
			case 1:	Nu[2] = GetRightU(t); break;
			case 2:	Mu[2] = 1; Nu[2] = h*GetRightU_x(t); break;
			case 3: Mu[2] = 1/(1 - h*rh); Nu[2] = -h*rh*GetRightTeta(t)/(1 - h*rh); break;
		}

		Progonka(N, A, B, C, F, Mu[1], Nu[1], Mu[2], Nu[2], U);
		for(i=0; i<=N; i++) Arr->get(j+1,i) = U[i];
	}
	
	delete []A;
	delete []B;
	delete []C;
	delete []F;
	delete []U;
}

void CEquation::SaveT(const char *fname, const double t)
{
	if(flag && t12->X1()<=t && t<=t12->X2())
	{
		int index = (int)( t/t12->H() );
		
		ofstream f(fname);
		for(int i=0; i<=x12->N(); i++) f<< x12->X(i) <<" "<< Arr->get(index,i) <<"\n";
		f.close();
	} else
	{
		cout << "CEquation::SaveT(): error\n";
		exit(EXIT_FAILURE);
	}
}

void CEquation::SaveX(const char *fname, const double x)
{
	if(x12->X1()<=x && x<=x12->X2())
	{
		int index = (int)( x/x12->H() );

		ofstream f(fname);
		for(int i=0; i<=t12->N(); i++) f<< t12->X(i) <<" "<< Arr->get(i,index) <<"\n";
		f.close();
	} else
	{
		cout << "CEquation::SaveX(): error\n";
		exit(EXIT_FAILURE);
	}
}