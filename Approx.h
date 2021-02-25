
/*
-------------------------------------------------------------------------------------------
	Цель:		Аппроксимация функции, заданной таблично, с заданным среднеквадратичым
				отклонением при помощи произвольного набора функций.

	Описание:	Конструктор: в М i-ый столбец содержит хi, yi, Roi соответственно; F(к,х) - значение к-й функции в х.

				Save: сохранение аппроксимации на произвольном интервале АВ.
				getDelta: среднеквадратичное отклонение.
				getF: значение аппроксиманты в произвольной точке.
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