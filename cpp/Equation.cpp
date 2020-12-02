
/*
------------------------------------------------------------------------------------------------------
	Файл:		Equation.cpp
	Версия:		1.11
	DLM:		17.08.2005
------------------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>
#include <math.h>

#include "Engine.h"
#include "Equation.h"
#include "TabFunc.h"

void SolveODE(CInterval AB, double U0, double (*fn)(double, double), char *fname)
{
	double
		x = AB.X(0),
		U = U0,
		h = AB.H(),
		k1,k2,k3,k4;

	ofstream f(fname); f<< x <<" "<< U <<"\n";

	for(int i=0; i<=AB.N(); i++)
	{
		x = AB.X(i);
		k1 = fn(x, U);
		k2 = fn(x + h/2., U + h*k1/2.);
		k3 = fn(x + h/2., U + h*k2/2.);
		k4 = fn(x + h, U + h*k3);

		f<< x <<" "<< (U += h/6.*(k1 + 2*(k2+k3) + k4)) <<"\n";
	}
	f.close();
}

CEqn::CEqn(double x1, double x2, double t2, int l, int r, double lH, double rH)
{
	assert(l>=1 && l<=3 && r>=1 && r<=3 && lH>0 && rH>0);

	D.x.ReBorn(x1, x2, 1);
	D.t.ReBorn(0, t2, 1);
	
	lbt = l; rbt = r; lh = lH; rh = rH;
}

void CEqn::sUt(char *fname, double t[], int size)
{
	for(int i=0; i<size; i++) assert(D.t.X1()<=t[i] && t[i]<=D.t.X2());

	ofstream f(fname);
	for(i=0; i<=D.x.N(); i++)
	{
		f<< D.x.X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr(D.t.i(t[j]),i);
		f<<"\n";
	}
	f.close();
}

void CEqn::sUx(char *fname, double x[], int size)
{
	for(int i=0; i<size; i++) assert(D.x.X1()<=x[i] && x[i]<=D.x.X2());

	ofstream f(fname);
	for(i=0; i<=D.t.N(); i++)
	{
		f<< D.t.X(i);
		for(int j=0; j<size; j++) f<<" "<< Arr(i,D.x.i(x[j]));
		f<<"\n";
	}
	f.close();
}

void CEqn::gUt(int it, CMatrix &m)
{
	int N = Arr.GetN();
	assert(0<=it && it<N);
	m.SetSize(2,N);

	for(int i=0; i<=N; i++)
	{
		m.x(i) = D.x.X(i);
		m.y(i) = Arr(it,i);
	}
}

void CEqn::gUx(int ix, CMatrix &m)
{
	int M = Arr.GetM();	
	assert(0<=ix && ix<M);
	m.SetSize(2,M);

	for(int i=0; i<=M; i++)
	{
		m.x(i) = D.t.X(i);
		m.y(i) = Arr(i,ix);
	}
}

void CEqn::gUt(double t, CMatrix &m) { gUt(D.t.i(t),m); }
void CEqn::gUx(double x, CMatrix &m) { gUx(D.x.i(x),m); }

void CEqn::Solve(double h, double tau)
{
	assert(h>0 && tau>0);						//установка шагов по пространственной и временной координатам
	D.x.ReBorn(D.x.X1(), D.x.X2(), h);
	D.t.ReBorn(D.t.X1(), D.t.X2(), tau);
	
	Arr.SetSize(D.t.N()+1, D.x.N()+1);							//Место для решения уравнения
	for(int i=0; i<=D.x.N(); i++) Arr(0,i) = gU0(D.x.X(i));		//задание начального значения
}

void CParEqn::Solve(double h, double tau)
{
	CEqn::Solve(h,tau);

	int N = D.x.N();
	CMatrix A(N),B(N),C(N),F(N),U(N+1);		//коэффициенты для метода прогонки
	double _2h2 = 2*h*h, _2h2_tau = _2h2/tau;

	for(int j=0; j<D.t.N(); j++)
	{
		for(int i=1; i<N; i++)
		{
			double
				_u = Arr(j,i-1),
				u = Arr(j,i),
				u_ = Arr(j,i+1),

				_x = D.x.X(i-1),
				x = D.x.X(i),
				x_ = D.x.X(i+1),
				t = D.t.X(j),

				Alpha =	(gK(x,t,u) + gK(x_,t,u_) + gV(x,t,u)*h)/2.,
				Beta =	(gK(x,t,u) + gK(_x,t,_u) - gV(x,t,u)*h)/2.,
				Gamma = _2h2_tau*gM(x,t,u);

			A(i) = Beta;
			B(i) = Alpha;
			C(i) = Alpha + Beta + Gamma;
			F(i) = u_*Alpha + _u*Beta - u*(Alpha + Beta - Gamma) + _2h2*gF(x,t,u);
		}
		
		int nj = j+1;
		double	Mu[3] = {0},
				Nu[3] = {0},
				t = D.t.X(nj);

		switch(lbt)
		{
			case 1:	Nu[1] = gLU(t);
				break;
			case 2:	Mu[1] = 1; Nu[1] = -h*gLdU_dx(t);
				break;
			case 3: Mu[1] = 1/(1 + h*lh); Nu[1] = h*lh*gLTeta(t)/(1 + h*lh);
				break;
		}
		switch(rbt)
		{
			case 1:	Nu[2] = gRU(t);
				break;
			case 2:	Mu[2] = 1; Nu[2] = h*gRdU_dx(t);
				break;
			case 3: Mu[2] = 1/(1 - h*rh); Nu[2] = -h*rh*gRTeta(t)/(1 - h*rh);
				break;
		}

		Progonka(N, A,B,C,F, Mu[1],Nu[1], Mu[2],Nu[2], U);
		for(i=0; i<=N; i++) Arr(nj,i) = U(i);
	}
}

void CHypEqn::Solve(double h, double tau)
{
	CEqn::Solve(h,tau);
	
	int N = D.x.N();
	CMatrix A(N),B(N),C(N),F(N),U(N+1);		//коэффициенты для метода прогонки

	double	_2h = 2*h,						//Для ускорения вычислений
			h2 = h*h,
			t_2 = tau/2.,
			h_2 = h/2.,
			h2_tau = h2/tau,
			_2h2_tau2 = 2*pow(h/tau,2);

	Arr(1,0) = gLU(tau);		//Задание граничных значений на первом слое
	Arr(1,N) = gRU(tau);

	//Вычисление значения функции на первом слое для запуска разностной схемы
	//
	for(int i=1; i<N; i++)
	{
		double
			_u = Arr(0,i-1),
			u = Arr(0,i),
			u_ = Arr(0,i+1),
			x = D.x.X(i);

		Arr(1,i) = u + tau*(gdU_dt(x) +	t_2/gM(x,0,u)*( 
							gK(x,0,u)/h2*(_u - 2*u + u_) + gV(x,0,u)/_2h*(u_ - _u) + gF(x,0,u) ));
	}

	//реализация разностной схемы
	//
	for(int j=0; j<=D.t.N()-2; j++)
	{
		for(i=1; i<N; i++)
		{
			double
				_u = Arr(j,i-1),
				u = Arr(j,i),
				u_ = Arr(j,i+1),

				x = D.x.X(i),
				t = D.t.X(j),

				Alpha =	gK(x,t,u) - gV(x,t,u)*h_2,
				Beta =	gK(x,t,u) + gV(x,t,u)*h_2,
				Gamma = h2_tau*gL(x,t,u),
				Delta = _2h2_tau2*gM(x,t,u);

			A(i) = Alpha;
			B(i) = Beta;
			C(i) = Alpha + Beta - Gamma + Delta;
			F(i) = _u*Alpha + u_*Beta - u*(Alpha + Beta + Gamma + Delta) + 2*(Arr(j+1,i)*Delta + gF(x,t,u)*h2);
		}

		int nj = j+2;
		double	Mu[3] = {0},
				Nu[3] = {0},
				t = D.t.X(nj);

		switch(lbt)
		{
			case 1:	Nu[1] = gLU(t); break;
			case 2:	Mu[1] = 1; Nu[1] = -h*gLdU_dx(t); break;
			case 3: Mu[1] = 1/(1 + h*lh); Nu[1] = h*lh*gLTeta(t)/(1 + h*lh); break;	//Надо бы ускорить (вынесением)
		}
		switch(rbt)
		{
			case 1:	Nu[2] = gRU(t); break;
			case 2:	Mu[2] = 1; Nu[2] = h*gRdU_dx(t); break;
			case 3: Mu[2] = 1/(1 - h*rh); Nu[2] = -h*rh*gRTeta(t)/(1 - h*rh); break;
		}

		Progonka(N, A,B,C,F, Mu[1],Nu[1], Mu[2],Nu[2], U);
		for(i=0; i<=N; i++) Arr(nj,i) = U(i);
	}
}