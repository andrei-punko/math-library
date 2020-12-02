
/*
------------------------------------------------------------------------------------------------
	����:		Equation.h
	������:		1.03
	DLM:		17.12.2003

	����:		������� ���������������� ���������

	��������:	1.������� ��� dU/dx = F(x,U) ������� �����-����� 4 �������.

				2.������� �� � ������� ����������� ��������������� ����. � ����������� ��������-
				�� ������� ����������� ������� (x1,x2,t2) � ���	��������� ������� (l,r). ����
				���� ������� 3-�� ����, ��������� (lH, rH) ���������� ���������� �����.
------------------------------------------------------------------------------------------------
*/

#ifndef EQUATION_H
#define EQUATION_H

#include "Engine.h"

//	������� ��� ������� �������
//
void SolveODU(CInterval AB, const double U0, double (*fn)(double, double), const char *fname);

//	������� �� ��������������� ����
//
class CEquation  
{
public:
	CEquation(const double x1, const double x2, const double t2,
		const int l = 1, const int r = 1, const double lH = 1, const double rH = 1);
	virtual ~CEquation();

	virtual double GetV(const double x, const double t, const double U) { return 0; };
	virtual double GetF(const double x, const double t, const double U) { return 0; };
	virtual double GetK(const double x, const double t, const double U) { return 0; };
	virtual double GetM(const double x, const double t, const double U) { return 0; };

	virtual double GetU0(const double x) { return 0; };
	
	virtual double GetLeftU(const double t) { return 0; };
	virtual double GetRightU(const double t) { return 0; };
	virtual double GetLeftU_x(const double t) { return 0; };
	virtual double GetRightU_x(const double t) { return 0; };
	virtual double GetLeftTeta(const double t) { return 0; };
	virtual double GetRightTeta(const double t) { return 0; };

	void Solve(const double h, const double tau);
	void SaveT(const char *fname, const double t);
	void SaveT(const char *fname, const double t[], const int size);

	void SaveX(const char *fname, const double x);
	void SaveX(const char *fname, const double x[], const int size);

	CMatrix<double>* CEquation::GetU_t(const double t);
	CMatrix<double>* CEquation::GetU_x(const double x);

private:
	CInterval *x12, *t12;
	CMatrix<double>	*Arr;
	bool flag;
	int lbt, rbt;
	double lh, rh;
};

#endif