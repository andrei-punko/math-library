
/*
-------------------------------------------------------------------------------------------
	Файл:		Approx.h
	Версия:		1.00
	DLM:		11.01.2004

	Цель:		Аппроксимация функции, заданной таблично, с заданным среднеквадратичым
				отклонением при помощи произвольного набора функций.

	Описание:	Конструктор: в массиве М i-ый столбец содержит хi, yi, Roi соответственно.
				F(к,х) - значение к-й функции в точке х.

				Save: сохранение аппроксимации на произвольном интервале АВ.
				getDelta: среднеквадратичное отклонение.
				getF: значение аппроксиманты в произвольной точке.
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