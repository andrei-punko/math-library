
/*
---------------------------------------------------------------------------------------------------
	Описание:	Численные значения констант и параметров веществ
---------------------------------------------------------------------------------------------------
*/

#ifndef CONSTANTS_H
#define	CONSTANTS_H

namespace ScCs				//	Scientific constants namespace
{
	const double
		PI = 3.14159,
		c = 3e+8,			//	Скорость света, м/с
		Mu0 = 1.25e-6,		//	Абсолютная магнитная проницаемость, Гн/м
		Eps0 = 8.85e-12,	//	Абсолютная диэлектрическая проницаемость, Ф/м
		N_A = 6.022e+23,	//	Число Авогадро, 1/моль
		k_B = 1.38e-23,		//	Постоянная Больцмана, Дж/К
		e_charge = 1.6e-19;	//	Заряд электрона, Кл
}

#endif