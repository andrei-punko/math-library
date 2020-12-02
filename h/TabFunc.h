
/*
------------------------------------------------------------------------------------------------
	����:		TabFunc.h
	������:		1.01
	DLM:		25.01.2004
	
	����:		�������, �������� ��������.
				����������� ����������� �� ����������������� (1 � 2 �����������), ��������������,
				������������, �������������.
				���������� ����� � �������� ������� ����� ���� ��������.
	
	��������:	�������������� - ������� �������� �� ����� � ���������� �����.
				����������������, �����������������, �������������� � �������������� ������� - ���
				������ ��������.
------------------------------------------------------------------------------------------------
*/

#ifndef TABFUNC_H
#define	TABFUNC_H

#include "Matrix.h"

class CTabFunc
{
public:
	CTabFunc(CMatrix &m);
	~CTabFunc();

	double D1(double x);
	double D2(double x);
	double Integrate();
	double Interpolate(double x);
	void Tabulate(char *fname);
	void Conversion(int size);

private:
	void CalcFactors();

	double *a, *b, *c, *d;
	int N;
	CMatrix M;
};

#endif