
#include <assert.h>
#include <math.h>

#include "Engine.h"

CInterval::CInterval(double X1, double X2, double H)
{
	ReBorn(X1, X2, H);
}

CInterval::CInterval(double X1, double X2, int N)
{
	ReBorn(X1, X2, N);
}

void CInterval::ReBorn(double X1, double X2, double H)
{
	assert(X1<X2 && H>0 && H<=X2-X1);
	x1 = X1; x2 = X2; h = H; n = (int)floor( (x2-x1)/h );	//если необходимо, n будет на 1 больше
}

void CInterval::ReBorn(double X1, double X2, int N)
{
	assert(X1<X2 && N>0);
	x1 = X1; x2 = X2; n = N; h = (x2-x1)/(double)n;
}

CInterval::CInterval(CInterval &interval)
{
	(*this) = interval;
}

CInterval &CInterval::operator=(CInterval &right)
{
	if(&right != this)
	{
		x1 = right.x1;
		x2 = right.x2;
		h = right.h;
		n = right.n;
	}
	return *this;
}