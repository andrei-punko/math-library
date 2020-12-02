
/*
----------------------------------------------------------------------------------------
	Файл:						Engine.cpp
	Версия:						1.06
	Дата последней модификации:	11.01.2004
----------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "Engine.h"
#include "Matrix.h"

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

double Simpson(CInterval AB, double (*getF)(double))
{
	if(0 != (AB.N())%2) AB.ReBorn(AB.X1(), AB.X2(), AB.N()+1);	//N должно быть четным

	double I = 0;
	for(int i=1; i<AB.N(); i+=2) I += getF( AB.X(i) );
	I *= 2;
	for(i=2; i<=AB.N()-2; i+=2) I += getF( AB.X(i) );
	I = AB.H()/3*(
		getF(AB.X1()) + I*2 + getF(AB.X2())
		);

	return I;
}

double Simpson(CMatrix &M, const double h)
{
	CInterval AB(M.x(0), M.x(M.GetN()-1), h);
	if(0 != AB.N()%2) AB.ReBorn(AB.X1(), AB.X2(), AB.N()+1);	//N должно быть четным

	CMatrix *nM = Conversion(M, AB.N()+1);
	double I = 0;
	for(int i=1; i<AB.N(); i+=2) I += nM->y(i);
	I *= 2;
	for(i=2; i<=AB.N()-2; i+=2) I += nM->y(i);
	I = h/3*( nM->y(0) + I*2 + nM->y(AB.N()) );
	delete nM;

	return I;
}

double Interpolation(CMatrix &M, const double x)
{									//Предполагается, что пары (X;Y) упорядочены по x
	int N = M.GetN()-1;
	
	assert(M.x(0)<=x && N>0 && x<=M.x(N));		//Проверка, принадлежности х интервалу интерполирования
	
	N++;
	double h1, h2,
	*A = new double[N],	*B = new double[N],
	*C = new double[N],	*F = new double[N],
	
	*a = new double[N],	*b = new double[N],
	*c = new double[N+1], *d = new double[N],
	*L = new double[N];
	assert(A!=0 && B!=0 && C!=0 && F!=0 && a!=0 && b!=0 && c!=0 && d!=0 && L!=0);
	
	N--;
	
	for(int i=1; i<N; i++)
	{
		h1 = A[i] = M.x(i) - M.x(i-1);
		h2 = B[i] = M.x(i+1) - M.x(i);

		C[i] = -2*(h1 + h2);
		F[i] = 3*( (M.y(i)-M.y(i-1))/h1 - (M.y(i+1)-M.y(i))/h2 );
	}
	Progonka(N, A, B, C, F, 0, 0, 0, 0, L);		//Методом прогонки решаем систему
	for(i=0; i<=N; i++) c[i+1] = L[i];
		
	delete []A; delete []B; delete []C; delete []F;	delete []L;

	for(i=1; i<=N; i++)
	{
		a[i] = M.y(i-1);
		h1 = M.x(i) - M.x(i-1);
		b[i] = (M.y(i) - M.y(i-1))/h1 - h1/3*(c[i+1] + 2*c[i]);
		d[i] = (c[i+1] - c[i])/(3*h1);
	}

	for(i=1; M.x(i)<x; i++);	//Поиск интервала [x_i-1; x_i], содержащего х

	h1 = x-M.x(i-1);
	h2 = a[i] + h1*(b[i] + h1*(c[i] + h1*d[i]));
	delete []a; delete []b; delete []c; delete []d;

	return h2;
}

CMatrix* Conversion(CMatrix &M, const int size)
{
	CMatrix *nM = new CMatrix(2, size);

	CInterval AB(nM->x(0) = M.x(0), nM->x(size-1) = M.x(M.GetN()-1), size-1);
	
	for(int i=0; i<size; i++)
	{
		nM->x(i) = AB.X(i);
		nM->y(i) = Interpolation(M, AB.X(i));
	}
	
	return nM;
}

void Save(CInterval &AB, double (*getF)(double), const char *fname)
{
	ofstream f(fname);
	for(int i=0; i<=AB.N(); i++) f<< AB.X(i) <<" "<< getF(AB.X(i)) <<"\n";
	f.close();
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

double CInterval::X(const int i)
{
	assert(i>=0 && i<=n);
	return x1 + i*h;
}

CMatrix* mSOLVE(CMatrix &A, CMatrix &B)
{
	int Size = A.GetM();
	assert(A.GetN() == Size && B.GetM() == Size);

	CMatrix *X = new CMatrix(Size);
	assert(X!=0);

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