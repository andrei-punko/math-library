
/*
-------------------------------------------------------------------------------------------
	����:		Approx.h
	������:		1.02
	DLM:		17.08.2005

	����:		������������� �������, �������� ��������, � �������� �����������������
				����������� ��� ������ ������������� ������ �������.

	��������:	�����������: � � i-�� ������� �������� �i, yi, Roi ��������������; F(�,�) - �������� �-� ������� � �.

				Save: ���������� ������������� �� ������������ ��������� ��.
				getDelta: ������������������ ����������.
				getF: �������� ������������� � ������������ �����.
-------------------------------------------------------------------------------------------
*/

#ifndef CAPPROX_H
#define CAPPROX_H

#include "Engine.h"

class CApprox
{
public:
	CApprox(CMatrix &m, double (*Fk)(int, double));
	~CApprox() {};
	
	void Approximate(int n=1);
	double getF(double x);
	double getDelta();
	void Save(CInterval &AB, char *fname);
	
private:
	CMatrix M,Ak;
	double (*getFk)(int, double);
	int n;
};

#endif