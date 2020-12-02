
/*
------------------------------------------------------------------------------------------------
	Файл:		TabFunc.cpp
	Версия:		1.01
	DLM:		25.01.2004
------------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>
#include "TabFunc.h"
#include "Engine.h"

CTabFunc::CTabFunc(CMatrix &m)
{
	assert(m.GetN()>1); M = m;	
	
	CalcFactors();
}

void CTabFunc::CalcFactors()
{
	N = M.GetN();
	double h1, h2,
	*A = new double[N],	*B = new double[N],
	*C = new double[N],	*F = new double[N],
	*L = new double[N];
	
	a = new double[N];	b = new double[N];
	c = new double[N+1]; d = new double[N];
	
	assert(A!=0 && B!=0 && C!=0 && F!=0 && a!=0 && b!=0 && c!=0 && d!=0 && L!=0);
	N--;
	
	for(int i=1; i<N; i++)
	{
		h1 = A[i] = M.x(i) - M.x(i-1);
		h2 = B[i] = M.x(i+1) - M.x(i);

		C[i] = -2.0*(h1 + h2);
		F[i] = 3.0*( (M.y(i)-M.y(i-1))/h1 - (M.y(i+1)-M.y(i))/h2 );
	}
	Progonka(N, A, B, C, F, 0, 0, 0, 0, L);		//Методом прогонки решаем систему
	for(i=0; i<=N; i++) c[i+1] = L[i];
		
	delete []A; delete []B; delete []C; delete []F;	delete []L;

	for(i=1; i<=N; i++)
	{
		a[i] = M.y(i-1);
		h1 = M.x(i) - M.x(i-1);
		b[i] = (M.y(i) - M.y(i-1))/h1 - h1/3.0*(c[i+1] + 2.0*c[i]);
		d[i] = (c[i+1] - c[i])/(3.0*h1);
	}
}

CTabFunc::~CTabFunc()
{
	delete []a;	delete []b;
	delete []c;	delete []d;
}

double CTabFunc::D1(double x)
{
	if((x<=M.x(0)) || (M.x(N)<=x)) return 0;
	
	for(int i=1; M.x(i)<x; i++);	//Поиск интервала [x_i-1; x_i], содержащего х

	double
		h1 = x-M.x(i-1),
		h2 = b[i] + h1*(2.0*c[i] + h1*3.0*d[i]);

	return h2;
}

double CTabFunc::D2(double x)
{
	if((x<=M.x(0)) || (M.x(N)<=x)) return 0;
	
	for(int i=1; M.x(i)<x; i++);	//Поиск интервала [x_i-1; x_i], содержащего х

	return 2.0*c[i] + 6.0*d[i]*(x-M.x(i-1));
}

double CTabFunc::Integrate()
{
	if(N%2 != 0) Conversion(N+1);
	double
		I = 0,
		h1, h2;

	for(int i=1; i<=N-1; i+=2)
	{
		h1 = 0.5*( M.x(i) - M.x(i-1) );
		h2 = 0.5*( M.x(i+1) - M.x(i) );

		I += M.y(i-1)*(3*h1-h2) + 4*M.y(i)*(h1+h2) + M.y(i+1)*(3*h2-h1);
	}

	return I/3.0;
}

double CTabFunc::Interpolate(double x)
{									//Предполагается, что пары (X;Y) упорядочены по x
	if(x<=M.x(0)) return M.y(0);
	if(M.x(N)<=x) return M.y(N);
	
	for(int i=1; M.x(i)<x; i++);	//Поиск интервала [x_i-1; x_i], содержащего х

	double
		h1 = x-M.x(i-1),
		h2 = a[i] + h1*(b[i] + h1*(c[i] + h1*d[i]));

	return h2;
}

void CTabFunc::Tabulate(char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=N; i++) f<< M.x(i) <<" "<< M.y(i) <<"\n";
	f.close();
}

void CTabFunc::Conversion(int size)
{
	assert(size>0);
	CMatrix nM(2, size);
	CInterval AB(M.x(0), M.x(N), size);
	
	for(int i=0; i<size; i++)
	{
		nM.x(i) = AB.X(i);
		nM.y(i) = Interpolate(AB.X(i));
	}

	M = nM; CalcFactors();
}