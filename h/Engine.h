
/*
------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.08
	DLM:		18.02.2005

	Цель:		-
	Описание:

	1. Класс CInterval:интервал [x1; x2], шаг h на нем, число кусков N, на которые он разбивается шагом h.
	2. Класс CMatrix: одно- или двумерный массив (размерность = кол-ву переменных в конструкторе).
	3. Метод прогонки.
	4. Интегрирование методом Симпсона.
	5. Интерполяция кубическим сплайном.
	6. Изменение (интерполированием) количества точек, задающих функцию (интервал постоянен).
	7. Табулирование функции на заданном интервале и сохранение результатов в файл.
	8. Нахождение минимума функции на заданном интервале методом золотого сечения.
	9. Решение системы линейных уравнений методом Гаусса с выбором главного элемента.
	10.Вычисление определителя матрицы (приведением ее к треугольной).
------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

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
	double X(int i);
	double H() { return h; };
	int N() { return n; };
	int i(const double x);

private:
	double x1, x2, h;
	int n;
};

//	Массив
//
class CMatrix
{
public:
	CMatrix(const int m = 1, const int n = 1);
	~CMatrix() { delete [] pData; };
	double& get(const int m, const int n = 0);

	const CMatrix &operator=(const CMatrix &);

	double& x(const int i) { return get(0,i); };
	double& y(const int i) { return get(1,i); };
	
	int GetM() { return M; };
	int GetN() { return N; };

	void Save(const char *fname);
	void Clear();
	void SwapLines(int, int);		//Перестановка строк
	void SwapCols(int, int);		//Перестановка столбцов

private:
	double *pData;
	int M, N;
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
double Simpson(CInterval &AB, double (*getF)(double));

//	Поиск минимума
//
double Min(CInterval &AB, double(*f)(double));

//	Метод Гаусса
//
CMatrix* mSOLVE(CMatrix &A, CMatrix &B);

//	Вычисление определителя
//
double mDET(CMatrix &A);

//	Табулирование функции на заданном интервале и сохранение результатов в файл
//
void Save(CInterval &AB, double (*func)(double), char *fname);

//	То же, но func(int) возвращает значение по индексу
//
void Save(CInterval &AB, double (*func)(int), char *fname);

//	Табулирование параметрической функции
//
void Save(CInterval &t, double (*x)(double), double (*func)(double), char *fname);

#endif