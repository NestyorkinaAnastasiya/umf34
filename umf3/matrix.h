/*matrix.h*/
#pragma once
#include <algorithm>
#include "basis.h"

namespace matrix
{
	struct Matrix
	{
		//������� ����
		//��������� ������ �����(��������)
		vector <int> ig;
		//������ �������� ��������������� ��������� 
		vector <int> jg;
		//������������ �������� �������
		vector <double> di;
		//��������������� �������� ������� ������������ �������
		vector <double> ggl;
		//��������������� �������� �������� ������������ �������
		vector <double> ggu;

		//����������� �������
		int n;
		//C����� ���������
		vector<vector<int>> adjacencyList;

		//��������� �������� �������
		void CreatePortret(int slaeSize, Grid grid);
		//������������� ������� ����� ��������� ��������
		void Initialize(int size1);
		//��������� ������� �� ������
		void MultiplyAx(const vector <double> a, vector<double>& result);
		//��������� ����������������� ������� �� ������
		void MultiplyATx(vector<double> a, vector<double>& result);


	};
}
