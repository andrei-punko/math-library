
/*
---------------------------------------------------------------------------------------------------
	Файл:		Equation.h
	Версия:		1.08
	DLM:		16.03.2005

	Цель:		Решение ДУ

	Описание:	1.Решение ОДУ dU/dx = F(x,U) методом Рунге-Кутта 4 порядка
				  Решение системы ОДУ dy_k/dx = f_k(x,y_1,...,y_k)
				  [ Не забывать break в разных case'ах функции f_k ! ]

				2.Решение параболического ДУ:
				M(x,t,U)*dU/dt = d/dx[K(x,t,U)*dU/dx] + V(x,t,U)*dU/dx + F(x,t,U)

				3.Решение гиперболического ДУ:
				M(x,t,U)*d^2U/dt^2 = K(x,t,U)*d^2U/dx^2 + V(x,t,U)*dU/dx + L(x,t,U)*dU/dt + F(x,t,U)
				[ Решение уравнения с ненулевыми функциями V и L не проверялось ]

				В конструктор передаются границы области (x1,x2,t2) и род граничных условий (l,r).
				Если есть условия 3-го рода, параметры (lH, rH) передаются последними двумя.
---------------------------------------------------------------------------------------------------
*/

#ifndef EQUATION_H
#define EQUATION_H

#include "Engine.h"

//	Решение ОДУ первого порядка
//
void SolveODU(CInterval AB, double U0, double (*fn)(double, double), char *fname);

//	Решение системы ОДУ первого порядка
//
void SolveODU(CInterval AB, double *U0, double (*fn)(double, double *, int), int N, char *fname);

//	Абстрактный класс, инкапсулирующий данные для ДУ 2 порядка в ЧП
//
class CEqn
{
public:
	struct Area { CInterval *x,*t; } D;
	CMatrix	*Arr;

	CEqn(double x1, double x2, double t2, int l=1, int r=1, double lH=1., double rH=1.);
	~CEqn();

	virtual double gU0(int ix) { return 0; }
	virtual double gLU(int it) { return 0; }
	virtual double gRU(int it) { return 0; }
	virtual double gLdU_dx(int it) { return 0; }
	virtual double gRdU_dx(int it) { return 0; }
	virtual double gLTeta(int it) { return 0; }
	virtual double gRTeta(int it) { return 0; }

	virtual double gM(int ix, int it, double U) { return 1.; }
	virtual double gK(int ix, int it, double U) { return 1.; }
	virtual double gV(int ix, int it, double U) { return 0; }
	virtual double gF(int ix, int it, double U) { return 0; }

	virtual void Solve(double h, double tau);

	void sUt(char *fname, double t) { sUt(fname, &t, 1); }
	void sUt(char *fname, double t[], int size);
	void sUx(char *fname, double x) { sUx(fname, &x, 1); }
	void sUx(char *fname, double x[], int size);

	CMatrix* gUt(int it);	//Получение среза U(x) при заданном t = t0
	CMatrix* gUx(int ix);	//Получение среза U(t) при заданном х = x0

protected:					//Доступ к элементам этой секции имеют только элементы класса и его потомки
	bool flag;
	int lbt, rbt;
	double lh, rh;

	virtual void plug() =0;	//Заглушка для невозможности создания элемента этого класса
};


//	Решение ДУ параболического типа
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

//	Решение ДУ гиперболического типа
//
class CHypEqn: public CEqn
{
public:
	CHypEqn(double x1, double x2, double t2, int l=1, int r=1, double lH=1, double rH=1): CEqn(x1,x2,t2,l,r,lH,rH) {};
	~CHypEqn() {};

	virtual double gdU_dt(int ix) { return 0; }
	virtual double gL(int ix, int it, double U) { return 0; };
	
	virtual void Solve(double h, double tau);

protected:
	virtual void plug() {};
};

#endif