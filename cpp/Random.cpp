
/*
---------------------------------------------------------------------------------------------------
	����:		Random.cpp
	������:		1.00
	DLM:		27.03.2005
---------------------------------------------------------------------------------------------------
*/

//	rand():			�������� ���� �� ����� "Numerical Recipes in C", �.279
//	rand_norm():	�������� ���� �� ����� ���� "��������� ����������������", ������ 3.4, �.4

#include <math.h>
#include <time.h>
#include <assert.h>

#include "Random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

long SEED;

void randomize()
{ 
	SEED = time(NULL);
}

double rand(double min, double max)
{
	assert(min<max);

	long k;
	double ans;

	SEED ^= MASK;
	k = SEED/MASK;
	SEED = IA*(SEED - k*IQ) - IR*k;
	if(SEED < 0) SEED += IM;
	ans = AM*SEED;
	SEED ^= MASK;

	return min + ans*(max - min);
}

double rand_norm(double Mu0, double D)
{
	assert(D>0);
	double v1,v2,s;
	
	do{
		v1 = 2*rand(0,1) - 1;
		v2 = 2*rand(0,1) - 1;
		s = pow(v1,2) + pow(v2,2);
	} while(s>=1);

	return Mu0 + D*v1*sqrt(-2*log(s)/s);
}