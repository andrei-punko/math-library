
/*
-------------------------------------------------------------------------------------------
	Файл:		Matrix.cpp
	Версия:		1.01
	DLM:		11.01.2004
-------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>

#include "Matrix.h"

CMatrix::CMatrix(const int m, const int n)
{
	assert(m>=1 && n>=1);
	pData = new double[(M = m)*(N = n)];
	assert(pData!=0);
}

double& CMatrix::get(const int m, const int n)
{
	assert(m<M && n<N && m>=0 && n>=0);
	return *(pData + m*N + n);
}

void CMatrix::Save(const char *fname)
{
	ofstream f(fname);
	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)	f<< get(i,j) <<" ";
		f << "\n";
	}
	f.close();
}

void CMatrix::SwapLines(int m1, int m2)
{
	assert(m1<M && m2<M && m1>=0 && m2>=0);
	
	double el;
	for(int i=0; i<N; i++)
	{
		el = *(pData + m1*N + i);
		*(pData + m1*N + i) = *(pData + m2*N + i);
		*(pData + m2*N + i) = el;
	}
}

void CMatrix::SwapCols(int n1, int n2)
{
	assert(n1<N && n2<N && n1>=0 && n2>=0);

	double el;
	for(int i=0; i<N; i++)
	{
		el = *(pData + i*N + n1);
		*(pData + i*N + n1) = *(pData + i*N + n2);
		*(pData + i*N + n2) = el;
	}
}

const CMatrix &CMatrix::operator=(const CMatrix &right)
{
	if(&right != this)
	{
		delete []pData;
		pData = new double[(M = right.M)*(N = right.N)];
		assert(pData!=0);

		for(int i=0; i<M; i++)
		for(int j=0; j<N; j++)
			*(pData + i*N + j) = *(right.pData + i*N + j);
	}
	return *this;
}

void CMatrix::Clear()
{
	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++)
		get(i,j) = 0;
}
