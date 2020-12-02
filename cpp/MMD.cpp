
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
	CMatrix F1(N,3), F2(N,3);	//������� ������ ���������� ��� ��� ���� ���������������� ��������� �����
	int i,k;					//��������
	double d = dt/(2*m_1), E_k = 0, P = 0;
	
	CalcF_i(F1);			//��������� ���������� ��� � ������������� ������� �������

	//���������� �������� �������� � ������� ������� � ��������� ������
	//
	ofstream f(fname);
	f<< "t Ek U W\n"<<
		"0 "<< (E_k = getE_k()) <<" "<< U <<" "<< (E_k + U) <<"\n";

	for(int step=1; step<=NSteps+NPSteps; step++)
	{
		CalcF_i(F2);		//��������� ���������� ��� �� ��������� ����

		for(k=0; k<3; k++)	//���������� ����� ���������� ��������� � ���������
		for(i=0; i<N; i++)	//���������� ������ ������
		{
			double dx = dt*(VXYZ->get(i,k) + d*F1.get(i,k));	//��������� �������� �������
			
			//� ������� ���� dt ������� ������ ��������� �� ����� L_2. � �������� ������ - ����������� ���������
			//
			if(dx > L_2) { cout<< "Step dt is large. Decrease it.\n"; exit(1); }

			XYZ->get(i,k) += dx;			//	����������� �������� ����������
			CheckBounds( XYZ->get(i,k) );	//	���������� �� ���� ������� �������� ����� ��������������� ������

			VXYZ->get(i,k) += d*(F1.get(i,k) + F2.get(i,k));	//��������� ����� �������� ���������
			F1.get(i,k) = F2.get(i,k);							//���������� ������ �� F2 � F1	//��������������!

			//	������ �������� (�.�. ����� ���������, ������ ��� ��� ���� 6 �����)
			//
			if(step > NSteps)
			{
				double x[2] = { XYZ->get(i,k)-L_2, XYZ->get(i,k)+L_2 };
				for(int y=0; y<2; y++) P -= x[y]*dV_dr( x[y] );
			}
		}

		if(step%10 == 0) cout<< step/10 <<endl;				//����������� ��������� �����

		f<< step*dt <<" "<< (E_k = getE_k()) <<" "<< U <<" "<< (E_k + U) <<"\n";	//���������� E_k,U,W
	}
	f.close();

	if(NPSteps != 0) P /= 6*NPSteps*L*L;

	return P;
}

//	������ ��� ������ ������� ��������� ���, ����������� �� ��� �� ������� ���� ��������� ������
//
void CLattice::CalcF_i(CMatrix &F)
{
	U = 0;		//	�������� ������������� ������� �������
	F.Clear();	//	���������� ���, ����������� �� ������ ������� �� ������� ���� ���������

	for(int i=0; i<N-1; i++)	//	�������� i-� �������
	for(int j=i+1; j<N; j++)	//	�������� � ��� � ���� ��� ������� � ����� �������� ��������
	{
		double
			d[3] = {0},			//	�������� ���������� ���
			R = 0;				//	� ���������� ����� ��������� i � j
		
		for(int k=0; k<3; k++)	//	���������� �������� ���������� (��) �� k-� ������������ ��� (��k)
		{
			d[k] = XYZ->get(i,k) - XYZ->get(j,k);	//	��k ����� i-� � j-� ���������
			CheckBounds( d[k] );					//	�� �������� ��k ������� ����������� �� ��k ����� i-� �������� � ������������� j-�

			R += pow( d[k],2 );	//	����������� ������� ���������� �� ������� ��k
		}
		R = sqrt(R);			//	��������� ���������� ����� ���������

		if(R <= L_2)			//	����� � ���� ���� ������ ������� � �������� ����� �������� L_2
		{
			U += V(R);				//	����� � ������������� ������� �� ��������������� ���� ������
			
			double df,				//	����� � k-� ���������� ����, ����������� �� i-� �������, �� ������� j-�
				dv_dr = dV_dr( R );
			for(k=0; k<3; k++)
			{						
				F.get(i,k) -= (df = d[k]*dv_dr);
				F.get(j,k) += df;	// 3-� ����� �������
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

	for(i=0; i<N-1; i++)		//	�������� i-� �������
	for(j=i+1; j<N; j++)		//	�������� � ��� � ���� ��� ������� � ����� �������� ��������
	{
		double
			d[3] = {0},			//	�������� ���������� ����
			R = 0;				//	� ���������� ����� ��������� i � j
		
		for(int k=0; k<3; k++)	//	���������� �������� ���������� (��) �� k-� ������������ ��� (��k)
		{
			d[k] = XYZ->get(i,k) - XYZ->get(j,k);	//	��k ����� i-� � j-� ���������
			CheckBounds( d[k] );					//	�� �������� ��k ������� ���. �� ��k ����� i-� �������� � ������������� j-�

			R += pow( d[k],2 );	//	����������� ������� ���������� �� ������� ��k
		}
		R = sqrt(R);			//	��������� ���������� ����� ���������
		M[AB.i(R)]++;			//	����������� �� ����������� �������, ��������������� ������� R
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

//	��������, �� �������� �� ������� �� ������� ����	[ ������ ���� ������� ���� �� �����; ��� ������ ������
//														���������� ����������� � Solve ���������� ������� dx <= L ]
//
inline void CLattice::CheckBounds(double &d)
{
//	if( fabs(d) > L_2 ) d -= SIGN(d)*L;		//������� ��������� �� ����� L �� dt
	d -= Round(d/L)*L;						//������� ��������� ����� ���������� �� dt
}

void CLattice::InitAmorph(double N_v)
{
	N = (int)(N_v*pow(L,3));
	Init();

	//	��������� ������� ��������� ������� � ����
	//
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) XYZ->get(i,k) = rand(-L_2,L_2);
}

void CLattice::InitVCC(double a)
{
	int N_line = floor(L/a),
		N__ = N_line-1;				//��� ��������� �����
	N = pow(N_line,3) + pow(N__,3);
	Init();

	//	��������� ������� �� ����� ���-�������
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
		D = sqrt(k_B*T/m_1),	//	��������� ��������� ���������
		V_mc[3] = {0};			//	���������� �������� ������ ���� (��)

	//	������ �������� ��������� ��������� ������ �� ��������� (�������, �.2, �.255), ��� ���� ��������� �������� ��
	//
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) V_mc[k] += (VXYZ->get(i,k) = rand_norm(0,D));

	//	��������� ��������, ����� �������� �� = 0
	//
	for(k=0; k<3; k++)
	{
		V_mc[k] /= (double)N;
		for(int i=0; i<N; i++) VXYZ->get(i,k) -= V_mc[k];
	}
}