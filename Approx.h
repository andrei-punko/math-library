
/*
-------------------------------------------------------------------------------------------
	����:		Approx.h
	������:		1.00
	DLM:		11.01.2004

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
	CApprox(CMatrix &m, Fk F, double Eps);
	~CApprox();
	
	double getF(double x, int n);
	double getF(double x) {	return getF(x, best.n); };
	
	double getDelta(int n);

	double getBdelta() { return best.delta; };
	int getBn() { return best.n; };

	void Save(CInterval &AB, const char *fname);

private:
	void Approximate(int n);
	
	CMatrix M, *Ak;
	Fk getFk;
	struct {
		int n;
		double delta;
	} best;
};

#endif