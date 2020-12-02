
/*
------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.08
	DLM:		18.02.2005

	����:		-
	��������:

	1. ����� CInterval:�������� [x1; x2], ��� h �� ���, ����� ������ N, �� ������� �� ����������� ����� h.
	2. ����� CMatrix: ����- ��� ��������� ������ (����������� = ���-�� ���������� � ������������).
	3. ����� ��������.
	4. �������������� ������� ��������.
	5. ������������ ���������� ��������.
	6. ��������� (�����������������) ���������� �����, �������� ������� (�������� ���������).
	7. ������������� ������� �� �������� ��������� � ���������� ����������� � ����.
	8. ���������� �������� ������� �� �������� ��������� ������� �������� �������.
	9. ������� ������� �������� ��������� ������� ������ � ������� �������� ��������.
	10.���������� ������������ ������� (����������� �� � �����������).
------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

//	��������
//
class CInterval
{
public:
	CInterval(const double X1, const double X2, const double H) { ReBorn(X1, X2, H); };
	CInterval(const double X1, const double X2, const int N) { ReBorn(X1, X2, N); };

	void ReBorn(const double X1, const double X2, const double H);
	void ReBorn(const double X1, const double X2, const int N);

	double X1() { return x1; };
	double X2() { return x2; };
	double X(int i);
	double H() { return h; };
	int N() { return n; };
	int i(const double x);

private:
	double x1, x2, h;
	int n;
};

//	������
//
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

//	��������
//
void Progonka(const int N,
			  const double *A, const double *B,
			  const double *C, const double *F,
			  const double m1, const double n1,
			  const double m2, const double n2,
			  double *Y);

//	��������������
//
double Simpson(CInterval &AB, double (*getF)(double));

//	����� ��������
//
double Min(CInterval &AB, double(*f)(double));

//	����� ������
//
CMatrix* mSOLVE(CMatrix &A, CMatrix &B);

//	���������� ������������
//
double mDET(CMatrix &A);

//	������������� ������� �� �������� ��������� � ���������� ����������� � ����
//
void Save(CInterval &AB, double (*func)(double), char *fname);

//	�� ��, �� func(int) ���������� �������� �� �������
//
void Save(CInterval &AB, double (*func)(int), char *fname);

//	������������� ��������������� �������
//
void Save(CInterval &t, double (*x)(double), double (*func)(double), char *fname);

#endif