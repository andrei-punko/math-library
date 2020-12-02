
/*
-------------------------------------------------------------------------------------------
	Файл:		Matrix.h
	Версия:		1.03
	DLM:		10.11.2003

	Цель:		Реализация одно- или двумерного массива переменных ЛЮБОГО типа.
				Размерность массива = количеству переменных в конструкторе

	Описание:	get возвращает ссылку на элемент массива, отслеживая выход за его
				границы. Результат get может быть L-value.
-------------------------------------------------------------------------------------------
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>		//Функция exit()

template <class T>
class CMatrix
{
public:
	CMatrix(const int m, const int n = 1);
	~CMatrix() { delete [] pData; };
	inline T& get(const int m, const int n = 0);

	T& x(const int i) { return get(0,i); };
	T& y(const int i) { return get(1,i); };
	
	int GetM() { return M; };
	int GetN() { return N; };

	void Save(const char *fname);
	void SwapLines(int, int);		//Перестановка строк
	void SwapCols(int, int);		//Перестановка столбцов

private:
	T *pData;
	int M, N;
};

template<class T>
CMatrix<T>::CMatrix(const int m, const int n)
{
	if(m<1 || n<1)
	{
		cout << "error matrix size\n";
		exit(EXIT_FAILURE);
	}
	pData = new T[(M = m)*(N = n)];
}

template<class T>
inline T& CMatrix<T>::get(const int m, const int n)
{
	if(m<M && n<N && m>=0 && n>=0) return *(pData + m*N + n);
	else
	{
		cout << "CMatrix::get(): error\n";
		exit(EXIT_FAILURE);
	}
}

template<class T>
void CMatrix<T>::Save(const char *fname)
{
	ofstream f(fname);
	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)	f<< get(i,j) <<" ";
		f << "\n";
	}
	f.close();
}

template<class T>
void CMatrix<T>::SwapLines(int m1, int m2)
{
	if(m1<M && m2<M && m1>=0 && m2>=0)
	{
		T el;

		for(int i=0; i<N; i++)
		{
			el = *(pData + m1*N + i);
			*(pData + m1*N + i) = *(pData + m2*N + i);
			*(pData + m2*N + i) = el;
		}
	} else
	{
		cout << "CMatrix::SwapLines(): error\n";
		exit(EXIT_FAILURE);
	}
}

template<class T>
void CMatrix<T>::SwapCols(int n1, int n2)
{
	if(n1<N && n2<N && n1>=0 && n2>=0)
	{
		T el;

		for(int i=0; i<N; i++)
		{
			el = *(pData + i*N + n1);
			*(pData + i*N + n1) = *(pData + i*N + n2);
			*(pData + i*N + n2) = el;
		}
	} else
	{
		cout << "CMatrix::SwapCols(): error\n";
		exit(EXIT_FAILURE);
	}
}

#endif
