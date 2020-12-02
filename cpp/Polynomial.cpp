
/*
-------------------------------------------------------------------------------------------
	Файл:		Polynomial.cpp
	Версия:		1.00
	DLM:		13.01.2004
-------------------------------------------------------------------------------------------
*/

#include "Polynomial.h"

//многочлен Лежандра
double polLejandr(int n, double x)
{
	switch(n)
	{
	case 0: return 1.0; break;
	case 1: return x; break;

	default: double N = (double)(--n);
		return (2.0*N+1)/(N+1)*x*polLejandr(n,x) - N/(N+1)*polLejandr(n-1,x);
	}
}

//многочлен Лагерра
double polLagerr(int n, double x)
{
	switch(n)
	{
	case 0: return 1.0; break;
	case 1: return x-1; break;

	default: double N = (double)(--n);
		return (x-2.0*N-1)*polLagerr(n,x) - N*N*polLagerr(n-1,x);
	}
}

//многочлен Эрмита
double polErmit(int n, double x)
{
	switch(n)
	{
	case 0: return 1.0; break;
	case 1: return 2.0*x; break;

	default: double N = (double)(--n);
		return 2.0*x*polErmit(n,x) - 2.0*N*polErmit(n-1,x);
	}
}

//многочлен Чебышева 1 рода
double polChebishev_1(int n, double x)
{
	switch(n)
	{
	case 0: return 1.0; break;
	case 1: return x; break;

	default: double N = (double)(--n);
		return 2.0*x*polChebishev_1(n,x) - polChebishev_1(n-1,x);
	}
}

//многочлен Чебышева 2 рода
double polChebishev_2(int n, double x)
{
	switch(n)
	{
	case 0: return 1.0; break;
	case 1: return 2.0*x; break;

	default: double N = (double)(--n);
		return 2.0*x*polChebishev_2(n,x) - polChebishev_2(n-1,x);
	}
}