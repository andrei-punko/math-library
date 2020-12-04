
/*
-----------------------------------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.12
	DLM:		17.08.2005
-----------------------------------------------------------------------------------------------------------------------------
	Цель:		-
	Описание:

	1.	Класс CInterval:интервал [x1; x2], шаг h на нем, число кусков N, на которые он разбивается шагом.

	2.	Класс CMatrix: одно- или двумерный массив (размерность = кол-ву переменных в конструкторе).

	3.	Метод прогонки.

	4.	Интегрирование методом Симпсона (интервал интегрирования может содержать нечетное количество частей).

	5.	Табулирование функции на заданном интервале и сохранение результатов в файл. 3 варианта реализации.

	6.	Нахождение минимума (максимума) функции на заданном интервале методом золотого сечения.

	7.	Решение системы линейных уравнений методом Гаусса с выбором главного элемента.

	8.	Нахождение обратной матрицы.

	9.	Вычисление определителя матрицы (приведением ее к треугольной).

	10.	Математическое округление.

	11.	Сохранение плотности распределения набора чисел с заданным шагом.

	12.	Сохранение плотности распределения функции на заданном интервале с заданным шагом.

	13.	Класс CLattice инкапсулирует параметры решетки частиц, расположенной в кубе:
			ребро куба и его половина			L, L_2(для ускорения вычислений)
			масса одной частицы					m_1
			количество частиц в кубе			N
			парный потенциал и его производная	V, dV_dr(для ускорения вычислений)
			координаты и скорости частиц		XYZ, VXYZ

		После создания объекта класса вызовом одной из ф-й Init...() выбирается тип решетки:
			InitAmorph(N_v)	аморфное в-во; частицы с концентрацей N_v равномерно размещаются в кубе случайным образом
			InitVCC(a)		ОЦК решетка с параметром а
			InitManual(m)	координаты каждого из атомов задаются "вручную"
		
		SaveKr: сохранение корелляционной ф. в fname; шаг гистограммы есть precision часть L: ширина столбика = precision*L

		SetTemperature(T)	установка скоростей частиц в соответствии с заданной температурой системы Т
		getE_k()			кинетическая энергия системы [не СРЕДНЕЕ. а суммирование кин. энергий частиц]
		getU()				потенциальная энергия системы

		По поводу реализации: для расчета сил, не теряя прежних значений, надо освобождать-выделять память, что долго
		!ВНИМАНИЕ: т.к. класс планировалось использовать только с ф. MMD, то данные в нем БЕЗРАЗМЕРНЫ; единицы выбраны
		следующим образом: единица пространственной координаты = периоду решетки, единица массы = массе атома.

  		CalcForces() расчитывает для каждой частицы компонент сил, действующих на нее со стороны всех остальных частиц.
		Учитывается взаимодействие только с частицами в пределах сферы r<L/2.
		После вызова CalcForces() обновляется значение lastU, которое и выдается ф-ей getU().
  		!НА БУДУЩЕЕ: предварительное табулирование V и dV_dr ускорит вычисления и избавит от передачи dV_dr в CLattice
	
	14.	Класс CUnSize реализует переход между знаковыми и беззнаковыми величинами. Масштаб задается набором единиц x,m,t.

	15.	Метод молекулярной динамики. Передаются длина шага, их количество, ссылка на подготовленный CLattice, флаг, функция
		для расчета средних. Флаг указывает на необходимость включения алгоритма, сохраняющего полную энергию. Функция для
		расчета средних определяется пользователем.
		!ВНИМАНИЕ: входные данные, потенциал и его производная ДОЛЖНЫ БЫТЬ БЕЗРАЗМЕРНЫМИ.
-----------------------------------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

#include <assert.h>
#include <math.h>

#include "Constants.h"

class CInterval
{
public:
	CInterval(double X1, double X2, double H);
	CInterval(double X1=0, double X2=1, int N=1);	//Конструктор по умолчанию
	CInterval(CInterval &);
	~CInterval() {};

	CInterval &operator=(CInterval &);

	void ReBorn(double X1, double X2, double H);
	void ReBorn(double X1, double X2, int N);

	inline double X1() { return x1; };
	inline double X2() { return x2; };
	inline double X(int i)
	{
		assert(0<=i && i<=n);
		return x1 + i*h;
	}

	inline double H() { return h; };
	inline int N() { return n; };
	inline int i(double x)
	{
		assert(x1<=x && x<=x2);
		return (int)( (x-x1)/h );
	}

private:
	double x1, x2, h;
	int n;
};

class CMatrix
{
public:
	CMatrix(int m=1, int n=1);
	CMatrix(CMatrix &);
	~CMatrix();

	CMatrix &operator=(CMatrix &);
	CMatrix &operator=(double d);
	
	CMatrix &operator-();
	bool operator==(CMatrix &);
	
	void SetSize(int m=1, int n=1);

	inline double& operator()(int m,int n=0)
	{
		assert(m>=0 && m<M && n>=0 && n<N);
		return *(pData + m*N + n);
	}

	inline double& x(int i) { return (*this)(0,i); };
	inline double& y(int i) { return (*this)(1,i); };
	
	inline int GetM() { return M; };
	inline int GetN() { return N; };

	double Max();
	double Min();
	
	void Save(char *fname, bool orig_view=true);
	void SwapLines(int, int);
	void SwapCols(int, int);

private:
	double *pData;
	int M, N;
};

class CUnSize
{
public:
	CUnSize(double a=1, double m=1, double tau=1): a(a),a2(a*a),m(m),tau(tau),tau2(tau*tau),alpha(m*pow(a/tau,2)),
		beta(m*pow(a/tau,2)/ScCs::k_B) {};
	~CUnSize() {};
	
	double X_s_to_u(double d=1) { return d/a; };
	double M_s_to_u(double d=1) { return d/m; };
	double Time_s_to_u(double d=1) { return d/tau; };
	double Vel_s_to_u(double d=1) { return tau/a*d; };
	double E_s_to_u(double d=1) { return d/alpha; };
	double T_s_to_u(double d=1) { return d/beta; };
	double S_s_to_u(double d=1) { return d/a2; };
	double F_s_to_u(double d=1) { return d*tau2/(m*a); };
	double P_s_to_u(double d=1) { return d*a*tau2/m; };
	
	double X_u_to_s(double d=1) { return a*d; };
	double M_u_to_s(double d=1) { return m*d; };
	double Time_u_to_s(double d=1) { return tau*d; };
	double Vel_u_to_s(double d=1) { return a/tau*d; };
	double E_u_to_s(double d=1) { return alpha*d; };
	double T_u_to_s(double d=1) { return beta*d; };
	double S_u_to_s(double d=1) { return a*a*d; };
	double F_u_to_s(double d=1) { return m*a/tau2*d; };
	double P_u_to_s(double d=1) { return m/(a*tau2)*d; };
	
private:
	double a,m,tau,alpha,beta, a2,tau2;
};

class CLattice
{
public:
	CLattice(double L, double (*V)(double), double (*dV_dr)(double)): L(L),L_2(L/2.),V(V),dV_dr(dV_dr),N(0) {};
	CLattice(CLattice &);
	~CLattice() {};
	
	CLattice &operator=(CLattice &);
	  
	void InitAmorph(double N_v);
	void InitVCC();
	void InitManual(CMatrix &XYZ);

	void SetTemperature(double T=0, double gamma=0.1);	//При выборе Т по умолчанию помнить: Т - уже безразмерная
	double getE_k();
	double getU();

	void SaveKr(char *fname, double precision=0.05);
	double CalcForces(CMatrix &F);
	
	double getRecommendedTau(double gamma=0.6);

	int N;
	double L,L_2,(*V)(double),(*dV_dr)(double);
	CMatrix XYZ,VXYZ;

private:
	double lastU;
};

//	Прогонка
//
void Progonka(int N, CMatrix &A,CMatrix &B,CMatrix &C,CMatrix &F, double m1,double n1,double m2,double n2, CMatrix &Y);

//	Интегрирование
//
double Simpson(CInterval &AB, double (*getF)(double));

//	Поиск минимума (максимума)
//
double FindMinMax(CInterval &AB, double(*f)(double), int find_max = -1);	//Плоховато

//	Метод Гаусса
//
void mSOLVE(CMatrix &A, CMatrix &B, CMatrix &X);

//	Нахождение обратной матрицы
//
void mINV(CMatrix &A, CMatrix &invA);

//	Вычисление определителя
//
double mDET(CMatrix &A);

//	Табулирование функции на заданном интервале и сохранение результатов в файл	//ИЗМЕНИТЬ НАЗВАНИЯ ФУНКЦИЙ !
//
void SaveFunc(CInterval &AB, double (*func)(double), char *fname);

//	То же, но func(int) возвращает значение по индексу
//
void SaveFunc(CInterval &AB, double (*func)(int), char *fname);

//	Табулирование параметрической функции
//
void SaveFunc(CInterval &t, double (*x)(double), double (*func)(double), char *fname);

//	Математическое округление
//
double Round(double d);

//	Сохранение плотности распределения (гистограммы) набора чисел с заданным шагом (шириной столбика)
//
void SaveCorrelationFunc(CMatrix &A, char *fname, double precision=0.05);

//	Сохранение плотности распределения (гистограммы) заданной функци на заданном интервале с заданным шагом (шириной столбика)
//
void SaveCorrelationFunc(CInterval &AB, double (*func)(double), char *fname, double precision=0.05);

//	))
//
void empty_func(CLattice &Lattice);

inline double sign(double d) { if(d<0) return -1; else return 1; }

//	ММД
//
void MMD(CLattice &Lattice, double dt, int StepsCount, bool sEnergy=false, void (*func)(CLattice &)=empty_func);

//	Если d не принадлежит [-L/2;L/2], оно изменяется на (+/-)L.
//	Если после этого d принадлежит интервалу, ф. возвращает true; или false в обратном случае.
//
inline bool CheckValue(double &d, double L)
{
	double L_2 = L/2.;
	if(fabs(d)>L_2) d -= L*sign(d);

	return fabs(d)<L_2;
}

//	То же для массива double величин
//
inline bool CheckValue(CMatrix &A, double L)
{
	for(int i=0; i<A.GetM(); i++)
	for(int j=0; j<A.GetN(); j++) if(!CheckValue(A(i,j),L)) return false;

	return true;
}

#endif