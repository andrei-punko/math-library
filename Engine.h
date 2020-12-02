
/*
------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.05
	DLM:		17.12.2003

	����:		-

	��������:
	1.����� CInterval ������������� �������� [x1; x2], ��� h �� ���, ����� ������ N,
		�� ������� �� ����������� ����� h.
	2.����� ��������.
	3.�������������� ������� ��������.
	4.������������ ���������� ��������.
	5.��������� (�����������������) ���������� �����, �������� ������� (�������� ���������).
	6.������������� ������� �� �������� ��������� � ���������� ����������� � ����.
	7.������� ������� �������� ��������� ������� ������ � ������� �������� ��������.
	8.���������� ������������ ������� (����������� �� � �����������).
------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

#include "Matrix.h"

//	��������
//
class CInterval
{
public:
	CInterval(const double X1, const double X2, const double H) { ReBorn(X1, X2, H); };
	CInterval(const double X1, const double X2, const int N) { ReBorn(X1, X2, N); };

	void ReBorn(const double X1, const double X2, const double H);
	void ReBorn(const double X1, const double X2, const int N);

	double X1() { return x1; };
	double X2() { return x2; };
	double X(const int i);
	double H() { return h; };
	int N() { return n; };

	int i(const double x);

private:
	double x1, x2, h;
	int n;
};

//	��������
//
void Progonka(const int N,
			  const double *A, const double *B,
			  const double *C, const double *F,
			  const double m1, const double n1,
			  const double m2, const double n2,
			  double *Y);

//	��������������
//
double Simpson(CMatrix &M, const double h);
double Simpson(CInterval AB, double (*getF)(double));

//	������������
//
double Interpolation(CMatrix &M, const double x);

//	��������� ���������� �����, �������� �������
//
CMatrix* Conversion(CMatrix &M, const int size);

//	�������������
//
void Save(CInterval &AB, double (*getF)(double), const char *fname);

//	����� ������
//
CMatrix* mSOLVE(CMatrix &A, CMatrix &B);

//	���������� ������������
//
double mDET(CMatrix &A);

#endif