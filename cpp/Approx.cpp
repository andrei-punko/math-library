
/*
-------------------------------------------------------------------------------------------
	Файл:		Approx.cpp
	Версия:		1.01
	DLM:		25.01.2004
-------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <math.h>
#include <fstream.h>
#include <iostream.h>

#include "Approx.h"
#include "Matrix.h"
#include "Engine.h"

CApprox::CApprox(CMatrix &m, Fk F):getFk(F)
{
	assert((m.GetM()==3) && (m.GetN()>1));
	M = m; Ak = 0;
}

void CApprox::Approximate(int _n)
{
	if(Ak!=0) delete Ak;

	n = _n;
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

double CApprox::getF(double x)
{
	double sum = 0;
	for(int k=0; k<=n; k++) sum += Ak->get(k)*getFk(k,x);

	return sum;
}

double CApprox::getDelta()
{
	double s1 = 0, s2 = 0;
	for(int i=0; i<M.GetN(); i++)
	{
		s1 += M.get(2,i)*pow(M.get(1,i) - getF(M.get(0,i)),2.0);
		s2 += M.get(2,i);
	}
	return pow(s1/s2, 0.5);
}

void CApprox::Save(CInterval &AB, char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< getF(AB.X(i)) <<"\n";
	f.close();
}