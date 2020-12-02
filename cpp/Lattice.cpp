
/*
------------------------------------------------------------------------------------------------
	File:		Lattice.cpp
	Version:	1.0
	DLM:		12.08.2005
------------------------------------------------------------------------------------------------
*/

#include <math.h>
#include <assert.h>
#include <iostream.h>	//�������!!!

#include "Engine.h"
#include "Random.h"
#include "Constants.h"

void CLattice::InitManual(CMatrix &XYZ)
{
	assert(XYZ.GetM()>=1 && XYZ.GetN()==3); (*this).XYZ = XYZ;

	SetTemperature();
}

void CLattice::InitAmorph(double N_v)
{
	N = (int)(N_v*pow(L,3));
	XYZ.SetSize(N,3);

	//	��������� ������� ��������� ������� � ����
	//
	my_randomize();

	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) XYZ(i,k) = my_rand(-L_2,L_2);

	SetTemperature();
}

void CLattice::InitVCC()
{
	int N_line = (int)floor(L),
		N__ = N_line-1;				//��� ��������� �����
	assert(N__>=1);

	N = (int)(pow(N_line,3) + pow(N__,3));
	XYZ.SetSize(N,3);

	//	��������� ������� �� ����� ���-�������
	//
	int m = 0;

	for(int x=0; x<N_line; x++)
	for(int y=0; y<N_line; y++)
	for(int z=0; z<N_line; z++)
	{
		XYZ(m,0) = -L_2 + x+0.5;
		XYZ(m,1) = -L_2 + y+0.5;
		XYZ(m,2) = -L_2 + z+0.5;
		m++;

		if(x!=N__ && y!=N__ && z!=N__)
		{
			XYZ(m,0) = -L_2 + x+1;
			XYZ(m,1) = -L_2 + y+1;
			XYZ(m,2) = -L_2 + z+1;
			m++;
		}
	}

	SetTemperature();
}

void CLattice::InitFCC()
{
	//�� ��������; ��������

	XYZ.SetSize(N=1,3);
	SetTemperature();
}

//
//	��������: gamma �� ������ ����� ������ ����� V!=0, �.� �� ������ ������ (int)(1/gamma); V0^2 = 3*T/gamma;
//	����������� - ���������� ��������
//
void CLattice::SetTemperature(double T, double gamma)
{
	assert(N>=1);
	VXYZ.SetSize(N,3);
	
	double V_mc[3] = {0},		//���������� �������� ������ ����
		V0 = sqrt(3*T/gamma);	//��������� �������� gamma ������ ����� ������
	
	my_randomize();
	int divider = (int)(1/gamma);

	//������ ��������� ��������, ��� ���� ��������� �������� ��
	//
	for(int i=0; i<N; i++)
	{
		double					//������� ����������� ����������
			phi = my_rand(0,ScCs::PI),
			teta = my_rand(0,ScCs::PI),
			V[3] = {sin(teta)*cos(phi), sin(teta)*sin(phi), cos(teta)};
		
		for(int k=0; k<3; k++) V_mc[k] += (VXYZ(i,k) = (i%divider ? 0 : V0*V[k]));
	}

	//	��������� ��������, ����� �������� �� = 0
	//
	for(int k=0; k<3; k++)
	{
		V_mc[k] /= (double)N;
		for(int i=0; i<N; i++) VXYZ(i,k) -= V_mc[k];
	}
}

CLattice &CLattice::operator=(CLattice &right)
{
	if(&right != this)
	{
		N = right.N;
		L = right.L;
		L_2 = right.L_2;
		V = right.V;
		dV_dr = right.dV_dr;
		XYZ = right.XYZ;
		VXYZ = right.VXYZ;
	}
	return *this;
}

CLattice::CLattice(CLattice &Lattice)
{
	(*this) = Lattice;
}

double CLattice::getE_k()
{
	double E_k = 0;
	
	for(int k=0; k<3; k++)
	for(int i=0; i<N; i++) E_k += pow(VXYZ(i,k),2);
		
	return E_k/2.;
}

void CLattice::SaveKr(char *fname, double precision)
{
	int Size = N*(N-1)>>1, count = 0;
	CMatrix Ri(Size);

	for(int i=0; i<N-1; i++)		//�������� i-� �������
	for(int j=i+1; j<N; j++)		//���������� ������� � �������� >i
	{
		double R=0, d[3]={0};		//�������� ���������� � ��� ����������

		for(int k=0; k<3; k++)		//���������� �������� ���������� (��) �� k-� ������������ ��� (��k)
		{
			d[k] = XYZ(i,k) - XYZ(j,k);	//��k ����� i-� � j-� ���������
			CheckValue(d[k],L);			//�������� ��k = ���. �� ��k ����� i-� �������� � ������������� j-�

			R += pow( d[k],2 );		//����������� ������� ���������� �� ������� ��������� ��������
		}
		
		Ri(count++) = sqrt(R);		//���������� ����� ���������
	}

	SaveCorrelationFunc(Ri,fname,precision);
}

double CLattice::getU()
{
	return lastU;
}

double CLattice::CalcForces(CMatrix &F)
{
	assert(F.GetM()==N && F.GetN()==3);

	F.Clear();						//F(i,k) - k-� ���������� ����, ����������� �� i-� ������� �� ������� ���������
	lastU = 0;						//������������� ������� �������

	for(int i=0; i<N-1; i++)		//�������� i-� �������
	for(int j=i+1; j<N; j++)		//���������� ������� � �������� j>i
	{
		double R=0, d[3]={0};		//�������� ���������� � ��� ����������

		for(int k=0; k<3; k++)		//���������� �������� ���������� (��) �� k-� ������������ ��� (��k)
		{
			d[k] = XYZ(i,k) - XYZ(j,k);	//��k ����� i-� � j-� ���������
			CheckValue(d[k],L);		//�������� ��k = ���. �� ��k ����� i-� �������� � ������������� j-�

			R += pow( d[k],2 );		//����������� ������� ���������� �� ������� ��������� ��������
		}
		R = sqrt(R);				//���������� ����� ���������

		if(R <= L_2)				//����� ���� ������ ������� � �������� ����� r<=L_2
		{
			lastU += V(R);			//����� � ������������� ������� �� ���� ������ i-j

			for(k=0; k<3; k++)
			{
				double df = d[k]/R*dV_dr(R);	//����� � F(i,k) �� ������� j-� �������
				F(i,k) -= df;
				F(j,k) += df;		//3-� ����� �������
			}
		}
	}

	return lastU;
}

double CLattice::getRecommendedTau(double gamma)
{
	assert(N>=1 && gamma>0 && gamma <1);

	double a = L/pow(N,1/3.),		//������ ��������� �������
			f = dV_dr(gamma*a);		//������ ������������ ���� [��������� �� gamma*�]

	return sqrt(2*(1-gamma)*a/fabs(f));
}