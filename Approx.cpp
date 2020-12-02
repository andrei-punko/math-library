
/*
-------------------------------------------------------------------------------------------
	Файл:		Approx.cpp
	Версия:		1.00
	DLM:		11.01.2004
-------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <math.h>
#include <fstream.h>
#include <iostream.h>

#include "CApprox.h"
#include "Matrix.h"
#include "Engine.h"

CApprox::CApprox(CMatrix &m, Fk F, double Eps):getFk(F)
{
	assert((m.GetM()==3) && (m.GetN()>1) && (Eps>0)); M = m;
	
	int n = 1;
	double d = 0;
	Approximate(n); best.n = n; best.delta = getDelta(n);
	
	do
	{
		Approximate(++n);
		d = getDelta(n);
		if(d<best.delta) { best.n = n; best.delta = d; }
	} while( (fabs(d/Eps-1)>1) && (n<M.GetN()) );

	if(n==M.GetN()) cout<< "\nCApprox: Approximation functions are not good.\n";
}

void CApprox::Approximate(int n)
{
	CMatrix A(n+1, n+1), B(n+1);
	A.Clear(); B.Clear();
	
	for(int p=0; p<=n; p++)	//цикл по строкам
	{
		for(int i=0; i<M.GetN(); i++)	//суммирование ряда
		{
			B.get(p) += M.get(2,i)*M.get(1,i)*getFk(p,M.get(0,i));

			for(int k=0; k<=n; k++)	//цикл по столбцам
				A.get(p, k) += M.get(2,i)*getFk(k, M.get(0,i))*getFk(p, M.get(0,i));
		}
	}
	Ak = mSOLVE(A, B);
}

CApprox::~CApprox()
{
	if(Ak!=0) delete Ak;
}

double CApprox::getF(double x, int n)
{
	double sum = 0;
	for(int k=0; k<=n; k++) sum += Ak->get(k)*getFk(k,x);

	return sum;
}

double CApprox::getDelta(int n)
{
	double s1 = 0, s2 = 0;
	for(int i=0; i<M.GetN(); i++)
	{
		s1 += M.get(2,i)*pow(M.get(1,i) - getF(M.get(0,i),n),2.0);
		s2 += M.get(2,i);
	}
	return pow(s1/s2, 0.5);
}

void CApprox::Save(CInterval &AB, const char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< getF(AB.X(i)) <<"\n";
	f.close();
}