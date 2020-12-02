
/*
------------------------------------------------------------------------------------------------
	Файл:		Equation.h
	Версия:		1.04
	DLM:		11.04.2004

	Цель:		Решение ДУ

	Описание:	1.Решение ОДУ dU/dx = F(x,U) методом Рунге-Кутта 4 порядка.

				2.Решение параболического ДУ:
				M(x,t,U)*dU/dt = d/dx[K(x,t,U)*dU/dx] + V(x,t,U)*dU/dx + F(x,t,U)
				
				В конструктор передаются границы области (x1,x2,t2) и род граничных условий (l,r).
				Если есть условия 3-го рода, параметры (lH, rH) передаются последними двумя.
------------------------------------------------------------------------------------------------
*/

#ifndef EQUATION_H
#define EQUATION_H

#include "Engine.h"

//	Решение ОДУ первого порядка
//
void SolveODU(CInterval AB, double U0, double (*fn)(double, double), char *fname);

//	Решение ДУ параболического типа
//
class CEquation  
{
public:
	struct { CInterval *x,*t; } D;
	CMatrix	*Arr;

	CEquation(double x1, double x2, double t2, int l=1, int r=1, double lH=1., double rH=1.);
	~CEquation();

	virtual double gV(int ix, int it, double U) { return 0; }
	virtual double gF(int ix, int it, double U) { return 0; }
	virtual double gK(int ix, int it, double U) { return 1.; }
	virtual double gM(int ix, int it, double U) { return 1.; }
	virtual double gU0(int ix) { return 0; }
	virtual double gLU(int it) { return 0; }
	virtual double gRU(int it) { return 0; }
	virtual double gLdU_dx(int it) { return 0; }
	virtual double gRdU_dx(int it) { return 0; }
	virtual double gLTeta(int it) { return 0; }
	virtual double gRTeta(int it) { return 0; }

	void Solve(double h, double tau);

	void sUt(char *fname, double t) { sUt(fname, &t, 1); }
	void sUt(char *fname, double t[], int size);
	void sUx(char *fname, double x) { sUx(fname, &x, 1); }
	void sUx(char *fname, double x[], int size);

	CMatrix* gUt(int it);
	CMatrix* gUx(int ix);

private:
	bool flag;
	int lbt, rbt;
	double lh, rh;
};

#endif