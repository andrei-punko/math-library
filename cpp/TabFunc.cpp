
/*
------------------------------------------------------------------------------------------------
	����:		TabFunc.cpp
	������:		1.03
	DLM:		17.08.2005
------------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>

#include "TabFunc.h"
#include "Engine.h"

CTabFunc::CTabFunc(CMatrix &m)
{
	assert(m.GetM()==2 && m.GetN()>=1);
	M = m;
	
	CalcFactors();
}

void CTabFunc::CalcFactors()
{
	N = M.GetN();
	
	double h1,h2;
	CMatrix A(N),B(N),C(N),F(N),L(N);
	
	a.SetSize(N); b.SetSize(N); c.SetSize(N+1); d.SetSize(N);
	N--;	//!?
	
	for(int i=1; i<N; i++)
	{
		h1 = A(i) = M.x(i) - M.x(i-1);
		h2 = B(i) = M.x(i+1) - M.x(i);

		C(i) = -2*(h1 + h2);
		F(i) = 3*((M.y(i) - M.y(i-1))/h1 - (M.y(i+1) - M.y(i))/h2);
	}
	Progonka(N, A, B, C, F, 0, 0, 0, 0, L);		//������� �������� ������ �������
	for(i=0; i<=N; i++) c(i+1) = L(i);

	for(i=1; i<=N; i++)
	{
		a(i) = M.y(i-1);
		h1 = M.x(i) - M.x(i-1);
		b(i) = (M.y(i) - M.y(i-1))/h1 - h1/3*(c(i+1) + 2*c(i));
		d(i) = (c(i+1) - c(i))/(3*h1);
	}
}

double CTabFunc::D1(double x)
{
	if(x<=M.x(0) || M.x(N)<=x) return 0;
	for(int i=1; M.x(i)<x; i++);		//����� ��������� [x_i-1; x_i], ����������� �

	double	h1 = x - M.x(i-1),
			h2 = b(i) + h1*(2*c(i) + h1*3*d(i));

	return h2;
}

double CTabFunc::D2(double x)
{
	if(x<=M.x(0) || M.x(N)<=x) return 0;
	for(int i=1; M.x(i)<x; i++);		//����� ��������� [x_i-1; x_i], ����������� �

	return 2*c(i) + 6*d(i)*(x - M.x(i-1));
}

double CTabFunc::D3(double x)
{
	if(x<=M.x(0) || M.x(N)<=x) return 0;
	for(int i=1; M.x(i)<x; i++);		//����� ��������� [x_i-1; x_i], ����������� �

	return 6*d(i);
}

double CTabFunc::Simpson()
{
	if(!(N%2)) Conversion(N+1);		//����������� �� ��������; �������� ������ �������� �� ������� ����� ������
	double h1,h2, I = 0;

	for(int i=1; i<N; i+=2)
	{
		h1 = 0.5*( M.x(i) - M.x(i-1) );
		h2 = 0.5*( M.x(i+1) - M.x(i) );

		I += M.y(i-1)*(3*h1-h2) + 4*M.y(i)*(h1+h2) + M.y(i+1)*(3*h2-h1);
	}

	return I/3.;
}

double CTabFunc::F(double x)
{
	if(x <= M.x(0)) return M.y(0);		//��������������, ��� ���� (X;Y) ����������� �� x
	if(M.x(N) <= x) return M.y(N);
	
	for(int i=1; M.x(i)<x; i++);		//����� ��������� [x_i-1; x_i], ����������� �

	double	h1 = x - M.x(i-1),
			h2 = a(i) + h1*(b(i) + h1*(c(i) + h1*d(i)));

	return h2;
}

void CTabFunc::SaveFunc(char *fname)
{
	M.Save(fname, false);
}

void CTabFunc::SaveFunc(CMatrix &m)
{
	m = M;
}

void CTabFunc::Conversion(int size)
{
	assert(size>=2);
	
	CMatrix nM(2, size);
	CInterval AB(M.x(0), M.x(N), size);
	
	for(int i=0; i<=size; i++)
	{
		nM.x(i) = AB.X(i);
		nM.y(i) = F(AB.X(i));
	}

	M = nM; CalcFactors();
}