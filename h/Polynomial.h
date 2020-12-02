
/*
-------------------------------------------------------------------------------------------
	����:		Polynomial.h
	������:		1.00
	DLM:		13.01.2004

	����:		���������� ��������, �������, ������, �������� 1 � 2 ����

	��������:	���������� ������������ ����� ������������ �������
-------------------------------------------------------------------------------------------
*/

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

//��������� ��������
double polLejandr(int n, double x);

//��������� �������
double polLagerr(int n, double x);

//��������� ������
double polErmit(int n, double x);

//��������� �������� 1 ����
double polChebishev_1(int n, double x);

//��������� �������� 2 ����
double polChebishev_2(int n, double x);

#endif