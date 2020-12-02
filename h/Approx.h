
/*
-------------------------------------------------------------------------------------------
	Файл:		Approx.h
	Версия:		1.01
	DLM:		25.01.2004

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