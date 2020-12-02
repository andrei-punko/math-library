
/*--------------------------------------------------------------------------------------
	Файл:						Engine.cpp
	Версия:						1.03
	Дата последней модификации:	10.11.2003
--------------------------------------------------------------------------------------*/

#include <stdlib.h>		//Функция exit()
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
	Alpha[1] = m1; Beta[1] = n1;

	for(int i=1; i<N; i++)
	{
		Alpha[i+1] = B[i]/(C[i]-A[i]*Alpha[i]);
		Beta[i+1] = (A[i]*Beta[i]+F[i])/(C[i]-A[i]*Alpha[i]);
	}

	Y[N] = (n2+m2*Beta[N])/(1-m2*Alpha[N]);
	for(i=N-1; i>=0; i--) Y[i] = Alpha[i+1]*Y[i+1]+Beta[i+1];

	delete Alpha;
	delete Beta;
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

double Simpson(CMatrix<double> &M, const double h)
{
	CInterval AB(M.x(0), M.x(M.GetN()-1), h);
	if(0 != AB.N()%2) AB.ReBorn(AB.X1(), AB.X2(), AB.N()+1);	//N должно быть четным

	CMatrix<double> *nM = Conversion(M, AB.N()+1);
	double I = 0;
	for(int i=1; i<AB.N(); i+=2) I += nM->y(i);
	I *= 2;
	for(i=2; i<=AB.N()-2; i+=2) I += nM->y(i);
	I = h/3*( nM->y(0) + I*2 + nM->y(AB.N()) );
	delete nM;

	return I;
}

double Interpolation(CMatrix<double> &M, const double x)
{									//Предполагается, что пары (X;Y) упорядочены по x
	int N = M.GetN()-1;
									//Проверка, принадлежности х интервалу интерполирования
	if(M.x(0)<=x && N>0 && x<=M.x(N))
	{
		N++;
		double h1, h2,
		*A = new double[N],	*B = new double[N],
		*C = new double[N],	*F = new double[N],
		
		*a = new double[N],	*b = new double[N],
		*c = new double[N+1], *d = new double[N],
		*L = new double[N];

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
		
		delete A; delete B; delete C; delete F;	delete L;

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
		delete a; delete b; delete c; delete d;

		return h2;
	} else
	{
		cout << "Interpolation(): error\n";
		exit(EXIT_FAILURE);
	}
}

CMatrix<double>* Conversion(CMatrix<double> &M, const int size)
{
	CMatrix<double> *nM = new CMatrix<double>(2, size);

	CInterval AB(nM->x(0) = M.x(0), nM->x(size-1) = M.x(M.GetN()-1), size-1);
	
	for(int i=0; i<size; i++)
	{
		nM->x(i) = AB.X(i);
		nM->y(i) = Interpolation(M, AB.X(i));
	}
	
	return nM;
}

void Save(CInterval AB, double (*getF)(double), const char *fname)
{
	ofstream f(fname);
	double
		x = AB.X1(),
		x2 = AB.X2(),
		h = AB.H();		

	while(x <= x2)
	{
		f<< x <<" "<< getF(x) <<"\n";
		x += h;
	}
	f.close();
}

void CInterval::ReBorn(const double X1, const double X2, const double H)
{
	if(X1<X2 && H>0 && H<=X2-X1)
	{
		x1 = X1; x2 = X2; h = H; n = (int)( (x2-x1)/h );
	} else
	{
		cout << "CInterval::SetH(): error\n";
		exit(EXIT_FAILURE);
	}
}

void CInterval::ReBorn(const double X1, const double X2, const int N)
{
	if(X1<X2 && N>0)
	{
		x1 = X1; x2 = X2; n = N; h = (x2-x1)/(double)n;
	} else
	{
		cout << "CInterval::SetH(): error\n";
		exit(EXIT_FAILURE);
	}
}

double CInterval::X(const int i)
{
	if(i>=0 && i<=n) return x1 + i*h;
	else
	{
		cout << "CInterval::X(): error\n";
		exit(EXIT_FAILURE);
	}
}

CMatrix<double>* mSOLVE(CMatrix<double> &A, CMatrix<double> &B)
{
	int Size = A.GetM();

	if(A.GetN() == Size && B.GetM() == Size)
	{
		CMatrix<double> *X = new CMatrix<double>(Size);

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
		for(int i=Size-1; i>=0; i--)
		{
			double d = B.get(i);
			for(int j=i+1; j<Size; j++) d -= A.get(i,j)*X->get(j);
			X->get(i) = d/A.get(i,i);
		}
		return X;
	} else
	{
		cout << "mSOLVE(): error\n";
		exit(EXIT_FAILURE);
	}
}

double mDET(CMatrix<double> &A)
{
	int Size = A.GetM();

	if(A.GetN() == Size)
	{
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
	} else
	{
		cout << "mDET(): error\n";
		exit(EXIT_FAILURE);
	}
}