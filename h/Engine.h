
/*
-----------------------------------------------------------------------------------------------------------------------------
	File:		Engine.h
	Version:	1.11
	DLM:		12.08.2005
-----------------------------------------------------------------------------------------------------------------------------
	����:		-
	��������:

	1. ����� CInterval:�������� [x1; x2], ��� h �� ���, ����� ������ N, �� ������� �� ����������� ����� h.
	2. ����� CMatrix: ����- ��� ��������� ������ (����������� = ���-�� ���������� � ������������).
	3. ����� ��������.
	4. �������������� ������� �������� (�������� �������������� ����� ��������� �������� ���������� ������).
	5. ������������� ������� �� �������� ��������� � ���������� ����������� � ����.
	6. ���������� �������� (���������) ������� �� �������� ��������� ������� �������� �������.
	7. ������� ������� �������� ��������� ������� ������ � ������� �������� ��������.
	8.���������� �������� �������.
	9.���������� ������������ ������� (����������� �� � �����������).
	10.�������������� ����������.
	11.���������� ��������� ������������� ������ ����� � �������� �����.
	12.���������� ��������� ������������� ������� �� �������� ��������� � �������� �����.
	13.�������� �������������� ������������ �������� ���������; ���������� � ��������� � ��������������� ������.
	14.
	����� CLattice ������������� ��������� ������� ������, ������������� � ����:
		����� ���� � ��� ��������			L, L_2(��� ��������� ����������)
		����� ����� �������					m_1
		���������� ������ � ����			N
		������ ��������� � ��� �����������	V, dV_dr(��� ��������� ����������)
		���������� � �������� ������		XYZ, VXYZ

	����� �������� ������� ������ ������� ����� �� �-� Init...() ���������� ��� �������:
		InitAmorph(N_v)	�������� �-��; ������� � ������������ N_v ���������� ����������� � ���� ��������� �������
		InitVCC(a)		��� ������� � ���������� �
		InitFCC(a)		��� ������� � ���������� � (�� ����������� !)
	
	SaveKr: ���������� �������������� �. � fname; ��� ����������� ���� precision ����� L: ������ �������� = precision*L

	SetTemperature(T)	��������� ��������� ������ � ������������ � �������� ������������ ������� �
	getE_k()			������������ ������� ������� [�� �������. � ������������ ���. ������� ������]
	getU()				������������� ������� �������

	�� ������ ����������: ��� ������� ���, �� ����� ������� ��������, ���� �����������-�������� ������, ��� �����

  	CalcForces() ����������� ��� ������ ������� ��������� ���, ����������� �� ��� �� ������� ���� ��������� ������.
	����������� �������������� ������ � ��������� � �������� ����� r<L/2.
	����� ������ CalcForces() ����������� �������� lastU, ������� � �������� �-�� getU().

  	!�� �������: ��������������� ������������� V � dV_dr ������� ���������� � ������� �� �������� dV_dr � CLattice
	
	15.����� CUnSize ��������� ������� ����� ��������� � ������������ ����������. ������� �������� ������� ������ x,m,t.
	16.����� ������������ ��������. ���������� ���, ���������� �����, ������ �� �������������� CLattice, ����, �������
	��� ������� �������. ���� ��������� �� ������������� ��������� ���������, ������������ ������ �������.
	!��������: ������� ������, ��������� � ��� ����������� ������ ���� �������������.
-----------------------------------------------------------------------------------------------------------------------------
*/

#ifndef ENGINE_H
#define ENGINE_H

#include <assert.h>
#include <math.h>

#include "Constants.h"

//	��������
//
class CInterval
{
public:
	CInterval(double X1, double X2, double H) { ReBorn(X1, X2, H); };
	CInterval(double X1, double X2, int N) { ReBorn(X1, X2, N); };
	CInterval(CInterval &);
	~CInterval() {};

	CInterval &operator=(CInterval &);

	void ReBorn(double X1, double X2, double H);
	void ReBorn(double X1, double X2, int N);

	inline double X1() { return x1; };
	inline double X2() { return x2; };
	inline double X(int i);
	inline double H() { return h; };
	inline int N() { return n; };
	inline int i(double x);

private:
	double x1, x2, h;
	int n;
};

//	������
//
class CMatrix
{
public:
	CMatrix(int m=1, int n=1);
	CMatrix(CMatrix &);
	~CMatrix();

	CMatrix &operator=(CMatrix &);
	CMatrix &operator=(double d);
	
	CMatrix &operator-();
	bool operator==(CMatrix &);
	
	void SetSize(int m=1, int n=1);

	inline double& operator()(int m,int n=0)
	{
		assert(m>=0 && m<M && n>=0 && n<N);
		return *(pData + m*N + n);
	}

	inline double& x(int i) { return (*this)(0,i); };
	inline double& y(int i) { return (*this)(1,i); };
	
	inline int GetM() { return M; };
	inline int GetN() { return N; };

	double Max();
	double Min();
	
	void Save(char *fname);
	void Clear();
	void SwapLines(int, int);
	void SwapCols(int, int);

private:
	double *pData;
	int M, N;
};

//	������� ������
//
class CLattice
{
public:
	CLattice(double L=1, double (*V)(double)=0, double (*dV_dr)(double)=0): L(L),L_2(L/2.),V(V),dV_dr(dV_dr),N(0) {};
	CLattice(CLattice &);
	~CLattice() {};
	
	CLattice &operator=(CLattice &);
	  
	void InitAmorph(double N_v);
	void InitVCC();
	void InitFCC();
	void InitManual(CMatrix &XYZ);

	void SetTemperature(double T=0, double gamma=0.1);	//��� ������ � �� ��������� �������: �- ��� ������������
	double getE_k();
	double getU();

	void SaveKr(char *fname, double precision=0.05);
	double CalcForces(CMatrix &F);
	
	double getRecommendedTau(double gamma=0.6);

	int N;
	double L,L_2,(*V)(double),(*dV_dr)(double);
	CMatrix XYZ,VXYZ;

private:
	double lastU;
};

//	��������
//
void Progonka(int N,
			  double *A, double *B, double *C, double *F,
			  double m1, double n1, double m2, double n2,
			  double *Y);

//	��������������
//
double Simpson(CInterval &AB, double (*getF)(double));

//	����� �������� (���������)
//
double FindMinMax(CInterval &AB, double(*f)(double), int find_max = -1);	//���������

//	����� ������
//
void mSOLVE(CMatrix &A, CMatrix &B, CMatrix &X);

//	���������� �������� �������
//
void mINV(CMatrix &A, CMatrix &invA);

//	���������� ������������
//
double mDET(CMatrix &A);

//	������������� ������� �� �������� ��������� � ���������� ����������� � ����	//�������� �������� ������� !
//
void Save(CInterval &AB, double (*func)(double), char *fname);

//	�� ��, �� func(int) ���������� �������� �� �������
//
void Save(CInterval &AB, double (*func)(int), char *fname);

//	������������� ��������������� �������
//
void Save(CInterval &t, double (*x)(double), double (*func)(double), char *fname);

//	�������������� ����������
//
double Round(double d);

//	���������� ��������� ������������� (�����������) ������ ����� � �������� ����� (������� ��������)
//
void SaveCorrelationFunc(CMatrix &A, char *fname, double precision=0.05);

//	���������� ��������� ������������� (�����������) �������� ������ �� �������� ��������� � �������� ����� (������� ��������)
//
void SaveCorrelationFunc(CInterval &AB, double (*func)(double), char *fname, double precision=0.05);

//	))
//
void empty_func(CLattice &Lattice);

inline double sign(double d) { if(d<0) return -1; else return 1; }

//	���
//
void MMD(CLattice &Lattice, double dt, int StepsCount, bool sEnergy=false, void (*func)(CLattice &)=empty_func);

//	�������� �������������� d [-L/2;L/2]. � ��������� ������ d ������������ � �������� ���������� �� nL (n - �����)
//
inline bool CheckValue(double &d, double L)
{
	double L_2 = L/2.;
	if(fabs(d)>L_2) d -= L*sign(d);

	return fabs(d)<L_2;
}

//	�� �� ��� �������
//
inline bool CheckValue(CMatrix &A, double L)
{
	for(int i=0; i<A.GetM(); i++)
	for(int j=0; j<A.GetN(); j++) if(!CheckValue(A(i,j),L)) return false;

	return true;
}

//	������� ����� ��������� � ������������ ����������
//
class CUnSize
{
public:
	CUnSize(double a, double m, double tau): a(a),a2(a*a),m(m),tau(tau),tau2(tau*tau),
		alpha(m*pow(a/tau,2)),beta(m*pow(a/tau,2)/ScCs::k_B) {};
	~CUnSize() {};

	double X_s_to_u(double d=1) { return d/a; };
	double M_s_to_u(double d=1) { return d/m; };
	double Time_s_to_u(double d=1) { return d/tau; };
	double Vel_s_to_u(double d=1) { return tau/a*d; };
	double E_s_to_u(double d=1) { return d/alpha; };
	double T_s_to_u(double d=1) { return d/beta; };
	double S_s_to_u(double d=1) { return d/a2; };
	double F_s_to_u(double d=1) { return d*tau2/(m*a); };
	double P_s_to_u(double d=1) { return d*a*tau2/m; };

	double X_u_to_s(double d=1) { return a*d; };
	double M_u_to_s(double d=1) { return m*d; };
	double Time_u_to_s(double d=1) { return tau*d; };
	double Vel_u_to_s(double d=1) { return a/tau*d; };
	double E_u_to_s(double d=1) { return alpha*d; };
	double T_u_to_s(double d=1) { return beta*d; };
	double S_u_to_s(double d=1) { return a*a*d; };
	double F_u_to_s(double d=1) { return m*a/tau2*d; };
	double P_u_to_s(double d=1) { return m/(a*tau2)*d; };

private:
	double a,m,tau,alpha,beta, a2,tau2;
};

#endif