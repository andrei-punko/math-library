
/*
-------------------------------------------------------------------------------------------
	Файл:		Approx.cpp
	Версия:		1.02
	DLM:		17.08.2005
-------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>
#include <math.h>

#include "Approx.h"
#include "Engine.h"

CApprox::CApprox(CMatrix &m, double (*Fk)(int, double))
{
	assert(m.GetM()==3 && m.GetN()>=2);

	getFk = Fk; M = m; Ak = 0;
}

void CApprox::Approximate(int n)
{
	assert(n>=1);

	CApprox::n = n;
	CMatrix A(n+1, n+1), B(n+1);
	A = 0; B = 0;					//Нельзя так: A = B = 0 (размеры матриц разные)
	
	for(int p=0; p<=n; p++)				//цикл по строкам
	{
		for(int i=0; i<M.GetN(); i++)	//суммирование ряда
		{
			B(p) += M(2,i)*M(1,i)*getFk(p,M(0,i));

			for(int k=0; k<=n; k++)		//цикл по столбцам
				A(p, k) += M(2,i)*getFk(k, M(0,i))*getFk(p, M(0,i));
		}
	}

	Ak.SetSize(n+1);
	mSOLVE(A,B,Ak);
}

double CApprox::getF(double x)
{
	double sum = 0;
	for(int k=0; k<=n; k++) sum += Ak(k)*getFk(k,x);

	return sum;
}

double CApprox::getDelta()
{
	double s1 = 0, s2 = 0;
	for(int i=0; i<M.GetN(); i++)
	{
		s1 += M(2,i)*pow(M(1,i) - getF(M(0,i)),2);
		s2 += M(2,i);
	}
	return pow(s1/s2, 0.5);
}

void CApprox::Save(CInterval &AB, char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< getF(AB.X(i)) <<"\n";
	f.close();
}