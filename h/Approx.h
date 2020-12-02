
/*
-------------------------------------------------------------------------------------------
	����:		Approx.h
	������:		1.01
	DLM:		25.01.2004

	����:		������������� �������, �������� ��������, � �������� �����������������
				����������� ��� ������ ������������� ������ �������.

	��������:	�����������: � ������� � i-�� ������� �������� �i, yi, Roi ��������������.
				F(�,�) - �������� �-� ������� � ����� �.

				Save: ���������� ������������� �� ������������ ��������� ��.
				getDelta: ������������������ ����������.
				getF: �������� ������������� � ������������ �����.
-------------------------------------------------------------------------------------------
*/

#ifndef CAPPROX_H
#define CAPPROX_H

#include "Matrix.h"
#include "Engine.h"

typedef double (*Fk)(int, double);

class CApprox
{
public:
	CApprox(CMatrix &m, Fk F);
	~CApprox();
	
	void Approximate(int _n);
	double getF(double x);
	double getDelta();
	void Save(CInterval &AB, char *fname);
	
private:
	CMatrix M, *Ak;
	Fk getFk;
	int n;
};

#endif