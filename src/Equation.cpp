#include <assert.h>
#include <fstream>
#include <math.h>

#include "Engine.h"
#include "Equation.h"
#include "TabFunc.h"

using namespace std;

void SolveODE(CInterval AB, double U0, double (*fn)(double, double), char *fname)
{
	double
		x = AB.X(0),
		U = U0,
		h = AB.H(),
		k1,k2,k3,k4;

	ofstream f(fname); f<< x <<" "<< U <<"\n";

	for(int i=0; i<=AB.N(); i++)
	{
		x = AB.X(i);
		k1 = fn(x, U);
		k2 = fn(x + h/2., U + h*k1/2.);
		k3 = fn(x + h/2., U + h*k2/2.);
		k4 = fn(x + h, U + h*k3);

		f<< x <<" "<< (U += h/6.*(k1 + 2*(k2+k3) + k4)) <<"\n";
	}
	f.close();
}