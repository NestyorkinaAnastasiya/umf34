#include "LU.h"
#include <math.h>
#include <iostream>
namespace lu
{
	LU::LU()
	{}

	LU::~LU()
	{}

	void LU::Convert(vector<int> _ig, vector<int> _jg, vector<double> _di, vector<double> _b, vector<double> _gl, vector<double> _gu, int size)
	{
		n = size;
		F.resize(n); D.resize(n); 
		x.resize(n); y.resize(n);
		
		ig.resize(n+1);
		for (int i = 0; i < n; i++)
		{
			F[i] = _b[i];
			D[i] = _di[i];
		}
		//��������� ����� ������ ig
		ig[0] = 0;
		for (int i = 1; i <= n; i++)
		{
			//���������� �� ������� ��������� � ������
			int k = _ig[i] - _ig[i - 1];
			//���� ����� ����
			if (k > 0) 
			{	// ��������� ������� ����� ��������� � ������, ��� � ����������
				int total_n = i - 1 - _jg[_ig[i-1]];
				ig[i] = ig[i - 1] + total_n;
			}
			else
				ig[i] = ig[i - 1];
		}
		U.resize(ig[n]);
		L.resize(ig[n]);

		// ��������� ����� gl � gu

		for (int i = 0; i < n; i++) 
		{	// ������ ������
			int start = ig[i];
			// �����
			int end = ig[i + 1];
			// ����� ������� �������
			int column = i - (end - start);
			int point = _ig[i];
			for (int j = start; j < end; j++, column++) 
			{
				if (column == _jg[point]) 
				{
					U[j] = _gu[point];
					L[j] = _gl[point];
					point++;
				}
				else U[j] = L[j] = 0;
			}
		}
	}

	void LU::Solve(vector<double> &u)
	{
		//���������� LU(sq) ����������
		Decomposition();

		//������� ������� Ly=b (������ ������) 
		LYF();

		//������� ������� Ux=y (�������� ������)
		UXY();
		u = x;
	}

	void LU::Decomposition()
	{	
		int i, j, i0, j0, iend, jend, kol_i, kol_j, num, first_i, first_j, ki, kj, dif;
		//����� ��� ��������� � L,U
		double suml, sumu, sumdg;

		for (i = 0; i < n; i++)
		{	//������ i ������ (�������)
			i0 = ig[i];
			//������ i+1 ������ (�������)
			iend = ig[i + 1];
			//���������� ��������� � i ������ (�������)
			kol_i = iend - i0;
			//����� ������� ���������� �������� i ������ (�������)
			first_i = i - kol_i;
			j = first_i;
			//��� �� ���� ��������� ��������� ������ (�������)
			for (num = i0, sumdg = 0; j<i; num++, j++)
			{
				j0 = ig[j];
				jend = ig[j + 1];
				kol_j = jend - j0;
				first_j = j - kol_j;
				ki = i0;
				kj = j0;
				//������� � ��������
				//�������� �� ������� ��� ����������� ���������
				dif = first_j - first_i;
				//ki - ����� ������� ����������� �������� � L,���� ��� �� ������
				if (dif > 0) ki += dif;
				//kj - ����� ������� ����������� �������� � U,���� ��� �� ������
				else kj -= dif;
				//��� num ����������� ��� ���������� ��������
				for (suml = 0, sumu = 0; ki<num; ki++, kj++) 
				{
					suml += L[ki] * U[kj];
					//��������� ����� ��� ������������� �������� �� U
					sumu += L[kj] * U[ki];
				}

				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[j];
				//���������� ������������ ��������
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}

	void LU::LYF()
	{	//// i0 - ����� ������ ������
		// iend - ����� ����� ������
		// ikol - ���������� ��������� � ������
		// beg - ����� ������� ������� ���������� �������� ������
		int i0, iend, ikol, beg; 
		double sum;

		for (int i = 0; i<n; i++)
		{
			i0 = ig[i]; iend = ig[i + 1];
			ikol = iend - i0;
			beg = i - ikol;
			sum = 0;
			for (int k = 0, j = beg; j<i; j++, k++)
				sum += y[j] * L[i0+k];

			y[i] = (F[i] - sum) / D[i];
		}

	}

	void LU::UXY()
	{	
		int i0, iend, ikol, beg;
		//������ �� �������� � �����
		for (int i = n - 1; i >= 0; i--)
		{
			x[i] += y[i];
			i0 = ig[i]; iend = ig[i + 1];
			ikol = iend - i0;
			beg = i - ikol;
			//��� �� ������� � �����
			for (int k = iend - 1, j = i - 1; j >= beg; j--, k--)
				x[j] -= x[i] * U[k];

		}

	}
}