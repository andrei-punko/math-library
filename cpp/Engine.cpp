
/*
----------------------------------------------------------------------------------------
	Файл:		Engine.cpp
	Версия:		1.11
	DLM:		14.03.2005
----------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <fstream.h>
#include <math.h>

#include "Engine.h"

#define L(x1,x2) (fabs(x1-x2))

void Progonka(const int N,
			  const double *A, const double *B,
			  const double *C, const double *F,
			  const double m1, const double n1,
			  const double m2, const double n2,
			  double *Y)
{
	double
		*Alpha = new double[N+1],
		*Beta = new double[N+1];
	assert(Alpha!=0 && Beta!=0);
	
	Alpha[1] = m1; Beta[1] = n1;

	for(int i=1; i<N; i++)
	{
		Alpha[i+1] = B[i]/(C[i]-A[i]*Alpha[i]);
		Beta[i+1] = (A[i]*Beta[i]+F[i])/(C[i]-A[i]*Alpha[i]);
	}

	Y[N] = (n2+m2*Beta[N])/(1-m2*Alpha[N]);
	for(i=N-1; i>=0; i--) Y[i] = Alpha[i+1]*Y[i+1]+Beta[i+1];

	delete []Alpha;
	delete []Beta;
}

double Simpson(CInterval &AB, double (*getF)(double))
{
	if(0 != (AB.N())%2) AB.ReBorn(AB.X1(), AB.X2(), AB.N()+1);	//N должно быть четным

	double I = 0;
	for(int i=1; i<AB.N(); i+=2) I += getF( AB.X(i) );
	I *= 2;
	for(i=2; i<=AB.N()-2; i+=2) I += getF( AB.X(i) );
	I = AB.H()/3.*(
		getF(AB.X1()) + I*2 + getF(AB.X2())
		);

	return I;
}

void CInterval::ReBorn(const double X1, const double X2, const double H)
{
	assert(X1<X2 && H>0 && H<=X2-X1);
	x1 = X1; x2 = X2; h = H; n = (int)( (x2-x1)/h );
}

void CInterval::ReBorn(const double X1, const double X2, const int N)
{
	assert(X1<X2 && N>0);
	x1 = X1; x2 = X2; n = N; h = (x2-x1)/(double)n;
}

double CInterval::X(int i)
{
	assert(i>=0 && i<=n);
	return x1 + i*h;
}

CMatrix* mSOLVE(CMatrix &a, CMatrix &b)
{
	int Size = a.GetM();
	assert(a.GetN() == Size && b.GetM() == Size);

	CMatrix
		*X = new CMatrix(Size),
		A(Size,Size), B(Size);
	assert(X!=0);

	A = a; B = b;	//Создаем копии матриц, т.к. над проводим действия

	//Приведение матрицы к треугольному виду
	for(int Curr=0; Curr<Size; Curr++)
	{
		int MaxNm = Curr;

		//Поиск максимального элемента в столбце Curr
		for(int i=Curr+1; i<Size; i++)
			if(A.get(i,Curr)>A.get(i,MaxNm)) MaxNm = i;

		//Перестановка строк
		if(MaxNm!=Curr)
		{
			A.SwapLines(Curr,MaxNm); B.SwapLines(Curr,MaxNm);
		}

		//Деление на коэффициент при первом члене
		B.get(Curr) = B.get(Curr)/A.get(Curr,Curr);
		for(i=Size-1; i>=Curr; i--)	A.get(Curr,i) = A.get(Curr,i)/A.get(Curr,Curr);			

		//Вычитание строк
		for(int line=Curr+1; line<Size; line++)
		{
			B.get(line) = B.get(line) - B.get(Curr)*A.get(line,Curr);
			for(i=Size-1; i>=Curr; i--)
				A.get(line,i) = A.get(line,i) - A.get(Curr,i)*A.get(line,Curr);
		}
	}

	//Нахождение решения системы по треугольной матрице
	double d;
	for(int i=Size-1; i>=0; i--)
	{
		d = B.get(i);
		for(int j=i+1; j<Size; j++) d -= A.get(i,j)*X->get(j);
		X->get(i) = d/A.get(i,i);
	}
		
	return X;
}

CMatrix* mINV(CMatrix &A)
{
	int Size = A.GetM();
	assert(Size == A.GetN());

	CMatrix
		E(Size),
		*invA = new CMatrix(Size,Size),
		*tM;
	
	for(int i=0; i<Size; i++)
	{
		E.Clear(); E.get(i) = 1;
		
		tM = mSOLVE(A,E);
		for(int j=0; j<Size; j++) invA->get(j,i) = tM->get(j);
		delete tM;
	}
	return invA;
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
		for(int i=Curr+1; i<Size; i++)
			if(A.get(i,Curr)>A.get(i,MaxNm)) MaxNm = i;

		//Перестановка строк
		if(MaxNm!=Curr)	{ A.SwapLines(Curr,MaxNm); NSw++; }
/*
		//Деление на коэффициент при первом члене
		for(i=Size-1; i>=Curr; i--)	A.get(Curr,i) = A.get(Curr,i)/A.get(Curr,Curr);			
*/
		//Вычитание строк
		for(int line=Curr+1; line<Size; line++)
		for(i=Size-1; i>=Curr; i--)
			A.get(line,i) = A.get(line,i) - A.get(Curr,i)/A.get(Curr,Curr)*A.get(line,Curr);
	}

	double DET = 1; if(NSw % 2 != 0) DET = -1;
	for(int i=0; i<Size; i++) DET *= A.get(i,i);

	return DET;
}

int CInterval::i(const double x)
{
	assert(x1<=x && x<=x2);
	return (int)( (x-x1)/h );
}

double Min(CInterval &AB, double(*f)(double))
{
	double x[4];
	x[0] = AB.X1(); x[2] = AB.X2();
	x[1] = x[0] + 2/(3+sqrt(5))*(x[2]-x[0]);

	while(x[2]-x[0]>AB.H())
	{
		x[3] = x[0]+x[2]-x[1];	
		int	m_i = 0, f_i = 0;
	
		//поиск минимальной из 4 точек
		for(int i=1; i<4; i++) if(f(x[i])<f(x[m_i])) m_i = i;

		//поиск наиболее удаленной от минимума точки
		for(i=1; i<4; i++) if(L(x[m_i],x[i]) > L(x[m_i],x[f_i])) f_i = i;

		//перенумеровываем точки, чтобы наиболее удаленная стала 4-й
		if(f_i!=3) x[f_i] = x[3];

		//упорядочение
		double d;

		for(m_i=0; m_i<=1; m_i++)
		{
			int min = m_i;
			for(i=1; i<=2; i++) if(x[i]<x[m_i]) min = i;

			d = x[min];
			x[min] = x[m_i];
			x[m_i] = d;
		}
	}
	return 0.5*(x[2]+x[0]);
}

void Save(CInterval &AB, double (*func)(int), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< func(i) <<"\n";
	f.close();
}

void Save(CInterval &AB, double (*func)(double), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< func(AB.X(i)) <<"\n";
	f.close();
}

void Save(CInterval &t, double(*x)(double), double (*func)(double), char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=t.N(); i++) f<< x( t.X(i) )<<" "<< func( t.X(i) ) <<"\n";
	f.close();
}

//---------------------------------------------------------------
//
CMatrix::CMatrix(const int m, const int n)
{
	assert(m>=1 && n>=1);
	pData = new double[(M = m)*(N = n)];
	assert(pData!=0);
}

double& CMatrix::get(const int m, const int n)
{
	assert(m<M && n<N && m>=0 && n>=0);
	return *(pData + m*N + n);
}

void CMatrix::Save(const char *fname)
{
	ofstream f(fname);
	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)	f<< get(i,j) <<" ";
		f << "\n";
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

const CMatrix &CMatrix::operator=(const CMatrix &right)
{
	if(&right != this)
	{
		delete []pData;
		pData = new double[(M = right.M)*(N = right.N)];
		assert(pData!=0);

		for(int i=0; i<M; i++)
		for(int j=0; j<N; j++)
			*(pData + i*N + j) = *(right.pData + i*N + j);
	}
	return *this;
}

void CMatrix::Clear()
{
	for(int i=0; i<M; i++)
	for(int j=0; j<N; j++) get(i,j) = 0;
}