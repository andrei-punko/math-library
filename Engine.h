
/*
------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.05
	DLM:		17.12.2003

	Цель:		-

	Описание:
	1.Класс CInterval инкапсулирует интервал [x1; x2], шаг h на нем, число кусков N,
		на которые он разбивается шагом h.
	2.Метод прогонки.
	3.Интегрирование методом Симпсона.
	4.Интерполяция кубическим сплайном.
	5.Изменение (интерполированием) количества точек, задающих функцию (интервал постоянен).
	6.Табулирование функции на заданном интервале и сохранение результатов в файл.
	7.Решение системы линейных уравнений методом Гаусса с выбором главного элемента.
	8.Вычисление определителя матрицы (приведением ее к треугольной).
------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

#include "Matrix.h"

//	Интервал
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

//	Прогонка
//
void Progonka(const int N,
			  const double *A, const double *B,
			  const double *C, const double *F,
			  const double m1, const double n1,
			  const double m2, const double n2,
			  double *Y);

//	Интегрирование
//
double Simpson(CMatrix &M, const double h);
double Simpson(CInterval AB, double (*getF)(double));

//	Интерполяция
//
double Interpolation(CMatrix &M, const double x);

//	Изменение количества точек, задающих функцию
//
CMatrix* Conversion(CMatrix &M, const int size);

//	Табулирование
//
void Save(CInterval &AB, double (*getF)(double), const char *fname);

//	Метод Гаусса
//
CMatrix* mSOLVE(CMatrix &A, CMatrix &B);

//	Вычисление определителя
//
double mDET(CMatrix &A);

#endif