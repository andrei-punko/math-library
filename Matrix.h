
/*
-------------------------------------------------------------------------------------------
	����:		Matrix.h
	������:		1.05
	DLM:		11.01.2004

	����:		���������� ����- ��� ���������� �������; ����������� ������� = ���������� ���������� �
				������������.

	��������:	get ���������� ������ �� ������� �������, ���������� ����� �� ���
				�������. ��������� get ����� ���� L-value.
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
	void SwapLines(int, int);		//������������ �����
	void SwapCols(int, int);		//������������ ��������

private:
	double *pData;
	int M, N;
};

#endif
