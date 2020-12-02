
/*
---------------------------------------------------------------------------------------------------
	Файл:		Random.cpp
	Версия:		1.02
	DLM:		09.08.2005
---------------------------------------------------------------------------------------------------
*/

//	my_rand():			Алгоритм взят из книги "Numerical Recipes in C", с.279
//	my_rand_norm():	Алгоритм взят из книги Кнут "Искусство программирования", раздел 3.4, с.4

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

void my_randomize()
{ 
	SEED = time(NULL);
}

double my_rand(double min, double max)
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

double my_rand_norm(double Mu0, double D)
{
	double v1,v2,s;
	
	do{
		v1 = 2*my_rand(0,1) - 1;
		v2 = 2*my_rand(0,1) - 1;
		s = pow(v1,2) + pow(v2,2);
	} while(s>=1);

	return Mu0 + D*v1*sqrt(-2*log(s)/s);
}