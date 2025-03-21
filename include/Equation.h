/*
---------------------------------------------------------------------------------------------------
	Цель:		Решение ДУ

	Описание:	1.Решение ОДУ dU/dx = F(x,U) методом Рунге-Кутта 4 порядка
---------------------------------------------------------------------------------------------------
*/

#ifndef EQUATION_H
#define EQUATION_H

#include "Engine.h"

//	Решение ОДУ первого порядка
//
void SolveODE(CInterval AB, double U0, double (*fn)(double, double), char *fname);

#endif