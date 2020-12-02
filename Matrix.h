
/*
-------------------------------------------------------------------------------------------
	Файл:		Matrix.h
	Версия:		1.05
	DLM:		11.01.2004

	Цель:		реализация одно- или двумерного массива; размерность массива = количеству переменных в
				конструкторе.

	Описание:	get возвращает ссылку на элемент массива, отслеживая выход за его
				границы. Результат get может быть L-value.
-------------------------------------------------------------------------------------------
*/

#ifndef MATRIX_H
#define MATRIX_H

class CMatrix
{
public:
	CMatrix(const int m = 1, const int n = 1);
	~CMatrix() { delete [] pData; };
	double& get(const int m, const int n = 0);

	const CMatrix &operator=(const CMatrix &);

	double& x(const int i) { return get(0,i); };
	double& y(const int i) { return get(1,i); };
	
	int GetM() { return M; };
	int GetN() { return N; };

	void Save(const char *fname);
	void Clear();
	void SwapLines(int, int);		//Перестановка строк
	void SwapCols(int, int);		//Перестановка столбцов

private:
	double *pData;
	int M, N;
};

#endif
