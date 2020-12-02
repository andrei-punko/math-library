
/*
---------------------------------------------------------------------------------------------------
	����:		Equation.h
	������:		1.09
	DLM:		17.08.2005

	����:		������� ��

	��������:	1.������� ��� dU/dx = F(x,U) ������� �����-����� 4 �������
				2.������� ��������������� �� M(x,t,U)*U'_t = [K(x,t,U)*U'_x]'_x + V(x,t,U)*U'_x + F(x,t,U)
				3.������� ���������������� �� M(x,t,U)*U''_tt = K(x,t,U)*U''_xx + V(x,t,U)*U'_x + L(x,t,U)*U'_t + F(x,t,U)
					! ��� V!=0, L!=0 �� �����������

				� ����������� ���������� ������� ������� (x1,x2,t2) � ��� ��������� ������� (l,r).
				���� ���� ������� 3-�� ����, ��������� (lH, rH) ���������� ���������� �����.
---------------------------------------------------------------------------------------------------
*/

#ifndef EQUATION_H
#define EQUATION_H

#include "Engine.h"

//	������� ��� ������� �������
//
void SolveODE(CInterval AB, double U0, double (*fn)(double, double), char *fname);

//	����������� �����, ��������������� ������ ��� �� � �� 2 �������
//
class CEqn
{
public:
	CEqn(double x1, double x2, double t2, int l=1, int r=1, double lH=1, double rH=1);
	~CEqn() {};

	virtual double gU0(double x) { return 0; }
	virtual double gLU(double t) { return 0; }
	virtual double gRU(double t) { return 0; }
	virtual double gLdU_dx(double t) { return 0; }
	virtual double gRdU_dx(double t) { return 0; }
	virtual double gLTeta(double t) { return 0; }
	virtual double gRTeta(double t) { return 0; }

	virtual double gM(double x, double t, double U) { return 1; }
	virtual double gK(double x, double t, double U) { return 1; }
	virtual double gV(double x, double t, double U) { return 0; }
	virtual double gF(double x, double t, double U) { return 0; }

	virtual void Solve(double h, double tau);

	void sUt(char *fname, double t) { sUt(fname, &t, 1); }
	void sUt(char *fname, double t[], int size);
	void sUx(char *fname, double x) { sUx(fname, &x, 1); }
	void sUx(char *fname, double x[], int size);

	void gUt(int it, CMatrix &m);
	void gUt(double t, CMatrix &);	//��������� ����� U(x) ��� �������� t
	void gUx(int ix, CMatrix &m);
	void gUx(double x, CMatrix &);	//��������� ����� U(t) ��� �������� x

protected:							//������ � ��������� ���� ������ ����� ������ �������� ������ � ��� �������
	int lbt, rbt;
	double lh, rh;
	struct Area{ CInterval x,t; } D;
	CMatrix	Arr;

	virtual void plug() =0;			//�������� ��� ������������� �������� �������� ����� ������
};

//	������� ��������������� ��
//
class CParEqn: public CEqn
{
public:
	CParEqn(double x1, double x2, double t2, int l=1, int r=1, double lH=1, double rH=1): CEqn(x1,x2,t2,l,r,lH,rH) {};
	~CParEqn() {};
	
	virtual void Solve(double h, double tau);

protected:
	virtual void plug() {};
};

//	������� ���������������� ��
//
class CHypEqn: public CEqn
{
public:
	CHypEqn(double x1, double x2, double t2, int l=1, int r=1, double lH=1, double rH=1): CEqn(x1,x2,t2,l,r,lH,rH) {};
	~CHypEqn() {};

	virtual double gdU_dt(double x) { return 0; }
	virtual double gL(double x, double t, double U) { return 0; };
	
	virtual void Solve(double h, double tau);

protected:
	virtual void plug() {};
};

#endif