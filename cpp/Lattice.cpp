
/*
------------------------------------------------------------------------------------------------
	File:		Lattice.cpp
	Version:	1.0
	DLM:		12.08.2005
------------------------------------------------------------------------------------------------
*/

#include <assert.h>
#include <math.h>

#include "../h/Engine.h"
#include "../h/Random.h"
#include "../h/Constants.h"

void CLattice::InitManual(CMatrix &XYZ)
{
	assert(XYZ.GetM()>=1 && XYZ.GetN()==3); (*this).XYZ = XYZ;

	SetTemperature();
}

void CLattice::InitAmorph(double N_v)
{
	N = (int)(N_v*pow(L,3));
	XYZ.SetSize(N,3);

	//	Случайным образом размещаем частицы в кубе
	//
	my_randomize();

	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) XYZ(i,k) = my_rand(-L_2,L_2);

	SetTemperature();
}

void CLattice::InitVCC()
{
	int N_line = (int)floor(L),
		N__ = N_line-1;				//Для ускорения счета
	assert(N__>=1);

	N = (int)(pow(N_line,3) + pow(N__,3));
	XYZ.SetSize(N,3);

	//	Размещаем частицы по узлам ОЦК-решетки
	//
	int m = 0;

	for(int x=0; x<N_line; x++)
	for(int y=0; y<N_line; y++)
	for(int z=0; z<N_line; z++)
	{
		XYZ(m,0) = -L_2 + x+0.5;
		XYZ(m,1) = -L_2 + y+0.5;
		XYZ(m,2) = -L_2 + z+0.5;
		m++;

		if(x!=N__ && y!=N__ && z!=N__)
		{
			XYZ(m,0) = -L_2 + x+1;
			XYZ(m,1) = -L_2 + y+1;
			XYZ(m,2) = -L_2 + z+1;
			m++;
		}
	}

	SetTemperature();
}

//
//	СКОРОСТИ: gamma от общего числа частиц имеет V!=0, т.е их номера кратны (int)(1/gamma); V0^2 = 3*T/gamma;
//	направления - выбираются случайно
//
void CLattice::SetTemperature(double T, double gamma)
{
	assert(N>=1);
	VXYZ.SetSize(N,3);
	
	double V_mc[3] = {0},		//Компоненты скорости центра масс
		V0 = sqrt(3*T/gamma);	//Начальная скорость gamma общего числа частиц
	
	my_randomize();
	int divider = (int)(1/gamma);

	//Задаем начальные скорости, при этом вычисляем скорость ЦМ
	//
	for(int i=0; i<N; i++)
	{
		double					//Введены сферические координаты
			phi = my_rand(0,ScCs::PI),
			teta = my_rand(0,ScCs::PI),
			V[3] = {sin(teta)*cos(phi), sin(teta)*sin(phi), cos(teta)};
		
		for(int k=0; k<3; k++) V_mc[k] += (VXYZ(i,k) = (i%divider ? 0 : V0*V[k]));
	}

	//	Нормируем скорости, чтобы скорость ЦМ = 0
	//
	for(int k=0; k<3; k++)
	{
		V_mc[k] /= (double)N;
		for(int i=0; i<N; i++) VXYZ(i,k) -= V_mc[k];
	}
}

CLattice &CLattice::operator=(CLattice &right)
{
	if(&right != this)
	{
		N = right.N;
		L = right.L;
		L_2 = right.L_2;
		V = right.V;
		dV_dr = right.dV_dr;
		XYZ = right.XYZ;
		VXYZ = right.VXYZ;
	}
	return *this;
}

CLattice::CLattice(CLattice &Lattice)
{
	(*this) = Lattice;
}

double CLattice::getE_k()
{
	double E_k = 0;
	
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) E_k += pow(VXYZ(i,k),2);
		
	return E_k/2.;
}

void CLattice::SaveKr(char *fname, double precision)
{
	int Size = N*(N-1)>>1, count = 0;
	CMatrix Ri(Size);

	for(int i=0; i<N-1; i++)		//Выбираем i-ю частицу
	for(int j=i+1; j<N; j++)		//Перебираем частицы с номерами >i
	{
		double R=0, d[3]={0};		//Обнуляем расстояние и его компоненты

		for(int k=0; k<3; k++)		//Перебираем проекции расстояний (ПР) на k-ю координатную ось (ПРk)
		{
			d[k] = XYZ(i,k) - XYZ(j,k);	//ПРk между i-й и j-й частицами
			CheckValue(d[k],L);			//Итоговое ПРk = мин. из ПРk между i-й частицей и отображениями j-й

			R += pow( d[k],2 );		//Увеличиваем квадрат расстояния на квадрат очередной проекции
		}
		
		Ri(count++) = sqrt(R);		//Расстояние между частицами
	}

	SaveCorrelationFunc(Ri,fname,precision);
}

double CLattice::getU()
{
	return lastU;
}

double CLattice::CalcForces(CMatrix &F)
{
	assert(F.GetM()==N && F.GetN()==3);

	F = 0;						//F(i,k) - k-я компонента силы, действующей на i-ю частицу со стороны остальных
	lastU = 0;						//Потенциальная энергия системы

	for(int i=0; i<N-1; i++)		//Выбираем i-ю частицу
	for(int j=i+1; j<N; j++)		//Перебираем частицы с номерами j>i
	{
		double R=0, d[3]={0};		//Обнуляем расстояние и его компоненты

		for(int k=0; k<3; k++)		//Перебираем проекции расстояний (ПР) на k-ю координатную ось (ПРk)
		{
			d[k] = XYZ(i,k) - XYZ(j,k);	//ПРk между i-й и j-й частицами
			CheckValue(d[k],L);		//Итоговое ПРk = мин. из ПРk между i-й частицей и отображениями j-й

			R += pow( d[k],2 );		//Увеличиваем квадрат расстояния на квадрат очередной проекции
		}
		R = sqrt(R);				//Расстояние между частицами

		if(R <= L_2)				//Вклад дают только частицы в пределах сферы r<=L_2
		{
			lastU += V(R);			//Вклад в потенциальную энергию от пары частиц i-j

			for(int k=0; k<3; k++)
			{
				double df = d[k]/R*dV_dr(R);	//Вклад в F(i,k) со стороны j-й частицы
				F(i,k) -= df;
				F(j,k) += df;		//3-й закон Ньютона
			}
		}
	}

	return lastU;
}

double CLattice::getRecommendedTau(double gamma)
{
	assert(N>=1 && gamma>0 && gamma <1);

	double a = L/pow(N,1/3.),		//Оценка параметра решетки
			f = dV_dr(gamma*a);		//Оценка максимальной силы [сближение на gamma*а]

	return sqrt(2*(1-gamma)*a/fabs(f));
}