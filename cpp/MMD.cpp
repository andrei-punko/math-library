
/*
------------------------------------------------------------------------------------------------
	File:		MMD.cpp
	Version:	1.01
	DLM:		28.06.2005
------------------------------------------------------------------------------------------------
*/

#include <fstream.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>	//exit()

#include "Engine.h"
#include "Random.h"
#include "Constants.h"
#include "MMD.h"

CLattice::~CLattice()
{
	if(!XYZ) delete XYZ;
	if(!VXYZ) delete VXYZ;
}

double CLattice::Solve(double dt, int NSteps, char *fname, int NPSteps)
{
	CMatrix F1(N,3), F2(N,3);	//Массивы хранят компоненты сил для двух последовательных временнЫх шагов
	int i,k;					//Счетчики
	double d = dt/(2*m_1), E_k = 0, P = 0;
	
	CalcF_i(F1);			//Вычисляем компоненты сил и потенциальную энергию системы

	//Записываем названия столбцов и энергии системы в начальный момент
	//
	ofstream f(fname);
	f<< "t Ek U W\n"<<
		"0 "<< (E_k = getE_k()) <<" "<< U <<" "<< (E_k + U) <<"\n";

	for(int step=1; step<=NSteps+NPSteps; step++)
	{
		CalcF_i(F2);		//Вычисляем компоненты сил на следующем шаге

		for(k=0; k<3; k++)	//Перебираем номер компоненты координат и скоростей
		for(i=0; i<N; i++)	//Перебираем номера частиц
		{
			double dx = dt*(VXYZ->get(i,k) + d*F1.get(i,k));	//Вычисляем смещение частицы
			
			//В течение шага dt частица должна пролетать не более L_2. В обратном случае - терминируем программу
			//
			if(dx > L_2) { cout<< "Step dt is large. Decrease it.\n"; exit(1); }

			XYZ->get(i,k) += dx;			//	Увеличиваем значение координаты
			CheckBounds( XYZ->get(i,k) );	//	Вылетевшие из куба частицы впускаем через противоположную стенку

			VXYZ->get(i,k) += d*(F1.get(i,k) + F2.get(i,k));	//Вычисляем новые значения скоростей
			F1.get(i,k) = F2.get(i,k);							//Сбрасываем данные из F2 в F1	//ОПТИМИЗИРОВАТЬ!

			//	Расчет давления (т.к. грани идентичны, делаем это для всех 6 сразу)
			//
			if(step > NSteps)
			{
				double x[2] = { XYZ->get(i,k)-L_2, XYZ->get(i,k)+L_2 };
				for(int y=0; y<2; y++) P -= x[y]*dV_dr( x[y] );
			}
		}

		if(step%10 == 0) cout<< step/10 <<endl;				//Отображение прогресса счета

		f<< step*dt <<" "<< (E_k = getE_k()) <<" "<< U <<" "<< (E_k + U) <<"\n";	//Записываем E_k,U,W
	}
	f.close();

	if(NPSteps != 0) P /= 6*NPSteps*L*L;

	return P;
}

//	Расчет для каждой частицы компонент сил, действующих на нее со стороны всех остальных частиц
//
void CLattice::CalcF_i(CMatrix &F)
{
	U = 0;		//	Обнуляем потенциальную энергию системы
	F.Clear();	//	Компоненты сил, действующих на каждую частицу со стороны всех остальных

	for(int i=0; i<N-1; i++)	//	Выбираем i-ю частицу
	for(int j=i+1; j<N; j++)	//	Выбираем к ней в пару все частицы с более высокими номерами
	{
		double
			d[3] = {0},			//	Обнуляем компоненты сил
			R = 0;				//	и расстояние между частицами i и j
		
		for(int k=0; k<3; k++)	//	Перебираем проекции расстояний (ПР) на k-ю координатную ось (ПРk)
		{
			d[k] = XYZ->get(i,k) - XYZ->get(j,k);	//	ПРk между i-й и j-й частицами
			CheckBounds( d[k] );					//	За итоговое ПРk считаем минимальное из ПРk между i-й частицей и отображениями j-й

			R += pow( d[k],2 );	//	Увеличиваем квадрат расстояния на квадрат ПРk
		}
		R = sqrt(R);			//	Вычисляем расстояние между частицами

		if(R <= L_2)			//	Вклад в силу дают только частицы в пределах сферы радиусом L_2
		{
			U += V(R);				//	Вклад в потенциальную энергию от рассматриваемой пары частиц
			
			double df,				//	Вклад в k-ю компоненту силы, действующей на i-ю частицу, со стороны j-й
				dv_dr = dV_dr( R );
			for(k=0; k<3; k++)
			{						
				F.get(i,k) -= (df = d[k]*dv_dr);
				F.get(j,k) += df;	// 3-й закон Ньютона
			}
		}
	}
}


void CLattice::SaveKr(double precision, char *fname)
{
	assert( precision>0 && precision<=1 );
	CInterval AB(0,L,L*precision);

	int *M = new int[AB.N()], i,j;
	for(i=0; i<AB.N(); i++) M[i] = 0;

	for(i=0; i<N-1; i++)		//	Выбираем i-ю частицу
	for(j=i+1; j<N; j++)		//	Выбираем к ней в пару все частицы с более высокими номерами
	{
		double
			d[3] = {0},			//	Обнуляем компоненты силы
			R = 0;				//	и расстояние между частицами i и j
		
		for(int k=0; k<3; k++)	//	Перебираем проекции расстояний (ПР) на k-ю координатную ось (ПРk)
		{
			d[k] = XYZ->get(i,k) - XYZ->get(j,k);	//	ПРk между i-й и j-й частицами
			CheckBounds( d[k] );					//	За итоговое ПРk считаем мин. из ПРk между i-й частицей и отображениями j-й

			R += pow( d[k],2 );	//	Увеличиваем квадрат расстояния на квадрат ПРk
		}
		R = sqrt(R);			//	Вычисляем расстояние между частицами
		M[AB.i(R)]++;			//	Увеличиваем на гистограмме столбик, соответствующий данному R
	}

	ofstream f(fname);
	for(i=0; i<AB.N(); i++)
	{
		double r = AB.H()*(i + 0.5);
		f<< r <<" "<< M[i]/(4*PI*pow(r,2)*AB.H()) <<"\n";
	}
	f.close();

	delete []M;
}

double CLattice::getE_k()
{
	double E_k = 0;
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) E_k += pow(VXYZ->get(i,k),2);

	return m_1*E_k/2.;
}

void CLattice::SavePositions(char *fname)
{
	XYZ->Save(fname);
}

void CLattice::SaveVelocities(char *fname)
{
	VXYZ->Save(fname);
}

void CLattice::Freeze()
{
	VXYZ->Clear();
}

//	Проверка, не вылетела ли частица за пределы куба	[ должна быть выбрана одна из строк; при выборе первой
//														необходимо отслеживать в Solve выполнение условия dx <= L ]
//
inline void CLattice::CheckBounds(double &d)
{
//	if( fabs(d) > L_2 ) d -= SIGN(d)*L;		//Частица пролетает не более L за dt
	d -= Round(d/L)*L;						//Частица пролетает любые расстояния за dt
}

void CLattice::InitAmorph(double N_v)
{
	N = (int)(N_v*pow(L,3));
	Init();

	//	Случайным образом размещаем частицы в кубе
	//
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) XYZ->get(i,k) = rand(-L_2,L_2);
}

void CLattice::InitVCC(double a)
{
	int N_line = floor(L/a),
		N__ = N_line-1;				//Для ускорения счета
	N = pow(N_line,3) + pow(N__,3);
	Init();

	//	Размещаем частицы по узлам ОЦК-решетки
	//
	int m = 0;

	for(int x=0; x<N_line; x++)
	for(int y=0; y<N_line; y++)
	for(int z=0; z<N_line; z++)
	{
		XYZ->get(m,0) = (a-L)/2. + a*x;
		XYZ->get(m,1) = (a-L)/2. + a*y;
		XYZ->get(m,2) = (a-L)/2. + a*z;
		m++;

		if(x!=N__ && y!=N__ && z!=N__)
		{
			XYZ->get(m,0) = -L/2. + a*(x+1);
			XYZ->get(m,1) = -L/2. + a*(y+1);
			XYZ->get(m,2) = -L/2. + a*(z+1);
			m++;
		}
	}
}

void CLattice::Init()
{
	XYZ = new CMatrix(N,3);
	VXYZ = new CMatrix(N,3);
	assert(XYZ && VXYZ);

	double
		D = sqrt(k_B*T/m_1),	//	Дисперсия гауссианы скоростей
		V_mc[3] = {0};			//	Компоненты скорости центра масс (ЦМ)

	//	Задаем проекции начальных скоростей частиц по Максвеллу (Сивухин, т.2, с.255), при этом вычисляем скорость ЦМ
	//
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) V_mc[k] += (VXYZ->get(i,k) = rand_norm(0,D));

	//	Нормируем скорости, чтобы скорость ЦМ = 0
	//
	for(k=0; k<3; k++)
	{
		V_mc[k] /= (double)N;
		for(int i=0; i<N; i++) VXYZ->get(i,k) -= V_mc[k];
	}
}