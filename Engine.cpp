
/*
----------------------------------------------------------------------------------------
	Файл:		Engine.cpp
	Версия:		1.14
	DLM:		17.08.2005
----------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <process.h>	//exit()

#include "Constants.h"
#include "Engine.h"

#define L(x1,x2)	(fabs(x1-x2))

using namespace std;

void Progonka(int N, CMatrix &A,CMatrix &B,CMatrix &C,CMatrix &F, double m1,double n1,double m2,double n2, CMatrix &Y)
{
	CMatrix Alpha(N+1), Beta(N+1);

	Alpha(1) = m1; Beta(1) = n1;
	for(int i=1; i<N; i++)
	{
		Alpha(i+1) = B(i)/(C(i) - A(i)*Alpha(i));
		Beta(i+1) = (A(i)*Beta(i) + F(i))/(C(i) - A(i)*Alpha(i));
	}

	Y(N) = (n2 + m2*Beta(N))/(1 - m2*Alpha(N));
	for(int i=N-1; i>=0; i--) Y(i) = Alpha(i+1)*Y(i+1) + Beta(i+1);
}

double Simpson(CInterval &AB, double (*getF)(double))
{
	if(AB.N() % 2) AB.ReBorn(AB.X1(), AB.X2(), AB.N()+1);	//N должно быть четным

	double I = 0;

	for(int i=1; i<AB.N(); i+=2) I += getF(AB.X(i));
	I *= 2;
	for(int i=2; i<=AB.N()-2; i+=2) I += getF(AB.X(i));
	I = AB.H()/3.*(
		getF(AB.X1()) + I*2 + getF(AB.X2())
		);

	return I;
}

void mSOLVE(CMatrix &a, CMatrix &b, CMatrix &X)
{
	int Size = a.GetM();
	assert(a.GetN() == Size && b.GetM() == Size);

	CMatrix	A(a), B(b);		//Создаем копии матриц, т.к. над ними проводим действия

	//Приведение матрицы к треугольному виду
	for(int Curr=0; Curr<Size; Curr++)
	{
		int MaxNm = Curr;

		//Поиск максимального элемента в столбце Curr
		for(int i=Curr+1; i<Size; i++)
			if(A(i,Curr)>A(i,MaxNm)) MaxNm = i;

		//Перестановка строк
		if(MaxNm!=Curr)
		{
			A.SwapLines(Curr,MaxNm); B.SwapLines(Curr,MaxNm);
		}

		//Деление на коэффициент при первом члене
		B(Curr) = B(Curr)/A(Curr,Curr);
		for(int i=Size-1; i>=Curr; i--)	A(Curr,i) = A(Curr,i)/A(Curr,Curr);

		//Вычитание строк
		for(int line=Curr+1; line<Size; line++)
		{
			B(line) = B(line) - B(Curr)*A(line,Curr);
			for(int i=Size-1; i>=Curr; i--)
				A(line,i) = A(line,i) - A(Curr,i)*A(line,Curr);
		}
	}

	//Нахождение решения системы по треугольной матрице
	double d;
	for(int i=Size-1; i>=0; i--)
	{
		d = B(i);
		for(int j=i+1; j<Size; j++) d -= A(i,j)*X(j);
		X(i) = d/A(i,i);
	}
}

void mINV(CMatrix &A, CMatrix &invA)
{
	assert(A.GetN() == A.GetM());

	int Size = A.GetM();
	invA.SetSize(Size,Size);

	CMatrix	E(Size),T(Size);
	for(int i=0; i<Size; i++)
	{
		E = 0; E(i) = 1;
		mSOLVE(A,E,T);

		for(int j=0; j<Size; j++) invA(j,i) = T(j);
	}
}

double mDET(CMatrix &A)
{
	int Size = A.GetM();
	assert(A.GetN() == Size);

	int NSw = 0; //Количество перестановок строк

	//Приведение матрицы к треугольному виду
	for(int Curr=0; Curr<Size; Curr++)
	{
		int MaxNm = Curr;

		//Поиск максимального элемента в столбце Curr
		for(int i=Curr+1; i<Size; i++) if(A(i,Curr)>A(i,MaxNm)) MaxNm = i;

		//Перестановка строк
		if(MaxNm!=Curr)	{ A.SwapLines(Curr,MaxNm); NSw++; }

		//Деление на коэффициент при первом члене
//		for(i=Size-1; i>=Curr; i--)	A(Curr,i) = A(Curr,i)/A(Curr,Curr);			

		//Вычитание строк
		for(int line=Curr+1; line<Size; line++)
		for(int i=Size-1; i>=Curr; i--)	A(line,i) = A(line,i) - A(Curr,i)/A(Curr,Curr)*A(line,Curr);
	}

	double DET = (NSw % 2) ? -1:1;
	for(int i=0; i<Size; i++) DET *= A(i,i);

	return DET;
}

double FindMinMax(CInterval &AB, double(*f)(double), int find_max)
{
	double x[4];
	x[0] = AB.X1(); x[2] = AB.X2();
	x[1] = x[0] + 2/(3+sqrt(5))*(x[2]-x[0]);

	while(x[2]-x[0]>AB.H())
	{
		x[3] = x[0]+x[2]-x[1];	
		int	m_i = 0, f_i = 0;
	
		//поиск минимальной из 4 точек
		for(int i=1; i<4; i++) if(find_max*(f(x[i])-f(x[m_i])) > 0) m_i = i;

		//поиск наиболее удаленной от минимума точки
		for(int i=1; i<4; i++) if(L(x[m_i],x[i]) > L(x[m_i],x[f_i])) f_i = i;

		//перенумеровываем точки, чтобы наиболее удаленная стала 4-й
		if(f_i!=3) x[f_i] = x[3];

		//упорядочение
		double d;

		for(m_i=0; m_i<=1; m_i++)
		{
			int min = m_i;
			for(int i=1; i<=2; i++) if(x[i]<x[m_i]) min = i;

			d = x[min];
			x[min] = x[m_i];
			x[m_i] = d;
		}
	}
	return 0.5*(x[2]+x[0]);
}

void SaveFunc(CInterval &AB, double (*func)(int), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< func(i) <<"\n";
	f.close();
}

void SaveFunc(CInterval &AB, double (*func)(double), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< func(AB.X(i)) <<"\n";
	f.close();
}

void SaveFunc(CInterval &t, double(*x)(double), double (*func)(double), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=t.N(); i++) f<< x( t.X(i) )<<" "<< func( t.X(i) ) <<"\n";
	f.close();
}

void SaveCorrelationFunc(CMatrix &A, char *fname, double precision)
{
	assert(0<precision && precision<=1 && A.GetM()>=1);
	
	double Min = A.Min(), Max = A.Max();
	CInterval AB(Min, Max, (Max-Min)*precision);
	CMatrix P(AB.N(),2);			//Частоты
	P = 0;
	
	for(int i=0; i<A.GetM(); i++)
	{
		double f_i = A(i);
		int num_of_f_i = AB.i(f_i);
		
		if(num_of_f_i != AB.N())	//CInterval включает обе границы диапазона 0..N; CMatrix - только 0
		{
			P(num_of_f_i,0) = A(i);
			P(num_of_f_i,1) += 1;
		}
	}
	
	//Нормировка
	//
	double max = P(0,1);
	for(int i=1; i<P.GetM(); i++) if(P(i,1) > max) max = P(i,1);
	for(int i=0; i<P.GetM(); i++) P(i,1) /= max;
	
	P.Save(fname);
}

void SaveCorrelationFunc(CInterval &AB, double (*func)(double), char *fname, double precision)
{
	CMatrix A(AB.N());
	for(int i=0; i<AB.N(); i++) A(i) = func(AB.X(i));
	
	SaveCorrelationFunc(A,fname,precision);
}

CMatrix::CMatrix(int m, int n)
{
	assert(m>=1 && n>=0);
	pData = new double[(M = m)*(N = n)];
	assert(pData);
}

void CMatrix::Save(char *fname, bool orig_view)
{
	ofstream f(fname);

	if(orig_view)
		for(int i=0; i<M; i++)
		{
			for(int j=0; j<N; j++)	f<< (*this)(i,j) <<" ";
			f<<"\n";
		}
	else
		for(int j=0; j<N; j++)
		{
			for(int i=0; i<M; i++)	f<< (*this)(i,j) <<" ";
			f<<"\n";
		}

	f.close();
}

void CMatrix::SwapLines(int m1, int m2)
{
	assert(m1<M && m2<M && m1>=0 && m2>=0);
	
	double el;
	for(int i=0; i<N; i++)
	{
		el = *(pData + m1*N + i);
		*(pData + m1*N + i) = *(pData + m2*N + i);
		*(pData + m2*N + i) = el;
	}
}

void CMatrix::SwapCols(int n1, int n2)
{
	assert(n1<N && n2<N && n1>=0 && n2>=0);

	double el;
	for(int i=0; i<N; i++)
	{
		el = *(pData + i*N + n1);
		*(pData + i*N + n1) = *(pData + i*N + n2);
		*(pData + i*N + n2) = el;
	}
}

CMatrix &CMatrix::operator=(CMatrix &right)
{
	if(&right != this)
	{
		SetSize(right.M,right.N);

		for(int i=0; i<M; i++)
		for(int j=0; j<N; j++) *(pData + i*N + j) = *(right.pData + i*N + j);
	}
	return *this;
}

CMatrix &CMatrix::operator=(double d)
{
	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) *(pData + i*N + j) = d;
		
	return *this;
}

CMatrix &CMatrix::operator-()
{
	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) *(pData + i*N + j) = -(*(pData + i*N + j));
		
	return *this;
}

bool CMatrix::operator==(CMatrix &right)
{
	assert(M==right.M && N==right.N);

	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) if(*(pData + i*N + j) != *(right.pData + i*N + j)) return false;
	
	return true;
}

double Round(double d)
{
	double r = ceil(d);		//Использование (int) недопустимо!
	if(fabs(d-r) > 0.5) return r += sign(d); else return r;
}

CMatrix::~CMatrix()
{
	if(!pData) delete [] pData;
}

void CMatrix::SetSize(int m, int n)		//!Память не освобождается, если L-массив УЖЕ имеет нужные размеры
{
	assert(m>=1 && n>=0);

	if(M!=m || N!=n)
	{
		if(!pData) delete [] pData;
		pData = new double[(M = m)*(N = n)];
		assert(pData);
	}
}

double CMatrix::Max()
{
	double Max = (*this)(0,0);
	
	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) if((*this)(i,j)>Max) Max = (*this)(i,j);

	return Max;
}

double CMatrix::Min()
{
	double Min = (*this)(0,0);

	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) if((*this)(i,j)<Min) Min = (*this)(i,j);

	return Min;
}

CMatrix::CMatrix(CMatrix &M)
{
	*this = M;
}

void empty_func(CLattice &Lattice) {}

void MMD(CLattice &Lattice, double dt, int StepsCount, bool sEnergy, void (*func)(CLattice &))
{
	int i,k, N = Lattice.N;
	CMatrix F1(N,3),F2(N,3);
	double t_2 = dt/2.;
	bool flag;
	
	Lattice.CalcForces(F1);
	
	double E0 = Lattice.getE_k() + Lattice.getU();	//Потенциальная энергия известна только после вызова CalcForces()

	for(int step=0; step<StepsCount; step++)
	{
		for(k=0; k<3; k++)
		for(i=0; i<N; i++) Lattice.XYZ(i,k) += dt*(Lattice.VXYZ(i,k) + t_2*F1(i,k));
		
		flag = CheckValue(Lattice.XYZ,Lattice.L);	//Вылетевшие частицы впускаем через противоположную стенку
		if(!flag)									//Проверка, не пролетает ли частица за dt больше L
		{
			cout<< "\nstep = "<< step <<"\n! dt is large !\n" <<endl; exit(1);
		}
		
		Lattice.CalcForces(F2);
			
		for(k=0; k<3; k++)
		for(i=0; i<N; i++) Lattice.VXYZ(i,k) += t_2*(F1(i,k) + F2(i,k));
		F1 = F2;

		if(sEnergy)
		{
			double H = 0;
			for(k=0; k<3; k++)
			for(i=0; i<N; i++) H += pow(Lattice.VXYZ(i,k),2) + pow(F2(i,k),2);

			double dE_H = (Lattice.getE_k()+Lattice.getU() - E0)/H;

			for(k=0; k<3; k++)
			for(i=0; i<N; i++)
			{
				Lattice.XYZ(i,k) += F2(i,k)*dE_H;
				Lattice.VXYZ(i,k) -= Lattice.VXYZ(i,k)*dE_H;
			}
		} //if
		
		func(Lattice);
	}
}