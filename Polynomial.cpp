
#include "Polynomial.h"

double polLejandr(int n, double x)
{
	switch(n)
	{
		case 0: return 1.0; break;
		case 1: return x; break;

		default: double N = (double)(--n); return (2.0*N+1)/(N+1)*x*polLejandr(n,x) - N/(N+1)*polLejandr(n-1,x);
	}
}

double polLagerr(int n, double x)
{
	switch(n)
	{
		case 0: return 1.0; break;
		case 1: return x-1; break;

		default: double N = (double)(--n); return (x-2.0*N-1)*polLagerr(n,x) - N*N*polLagerr(n-1,x);
	}
}

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

double polChebishev_1(int n, double x)
{
	switch(n)
	{
		case 0: return 1.0; break;
		case 1: return x; break;

		default: 
			return 2.0*x*polChebishev_1(n-1,x) - polChebishev_1(n-2,x);
	}
}

double polChebishev_2(int n, double x)
{
	switch(n)
	{
		case 0: return 1.0; break;
		case 1: return 2.0*x; break;

		default: 
			return 2.0*x*polChebishev_2(n-1,x) - polChebishev_2(n-2,x);
	}
}