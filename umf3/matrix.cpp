/*matrix.cpp*/
#include "matrix.h"

namespace matrix
{
	//��������� �������� �������
	void Matrix::CreatePortret(int slaeSize, Grid grid)
	{
		n = slaeSize * 2;
		adjacencyList.resize(slaeSize * 2);
		int node1, node2;
		bool flag = false;
		for (int i = 0; i < grid.elements.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < j; k++)
				{
					node1 = grid.elements[i].nodes[j];
					node2 = grid.elements[i].nodes[k];
					if (node1 < node2)
						swap(node1, node2);
					flag = false;
					for (auto it = adjacencyList[node1 * 2].begin(); it != adjacencyList[node1 * 2].end(); ++it)
					{
						if (*it == 2 * node2) {
							flag = true; break;
						}
					}
					if (!flag)
					{
						adjacencyList[node1 * 2].push_back(node2 * 2);
						adjacencyList[2 * node1 + 1].push_back(node2 * 2);
						adjacencyList[node1 * 2].push_back(node2 * 2 + 1);
						adjacencyList[2 * node1 + 1].push_back(node2 * 2 + 1);
					}

				}
			}
		}
		for (int i = 0; i < slaeSize; i++)
			adjacencyList[2 * i + 1].push_back(2 * i);
		//// ���������� ������ ��������� �� �����������
		//sort(adjacencyList.begin(), adjacencyList.end());

		/* ���������� */
		double tmp;
		int elem, size;
		for (int elem = 0; elem < 2 * slaeSize; elem++)
		{
			size = adjacencyList[elem].size();
			// i - ����� �������
			for (int i = 0; i < size - 1; ++i) 
			{	// ���������� ���� �������
				for (int j = 0; j < size - 1; ++j) 
				{
					if (adjacencyList[elem][j + 1] < adjacencyList[elem][j])
					{
						tmp = adjacencyList[elem][j + 1];
						adjacencyList[elem][j + 1] = adjacencyList[elem][j];
						adjacencyList[elem][j] = tmp;
					}
				}
			}
		}
		
		/* ������������� ������� */
		Initialize(slaeSize);
}

	//������������� ������� ����� ��������� ��������
	void Matrix::Initialize(int slaeSize)
	{
		int size;
		ig.resize(2 * slaeSize + 1);
		di.resize(2 * slaeSize);

		ig[0] = 0;
		ig[1] = 0;
		for (int i = 1; i < 2 * slaeSize; i++)
		{
			size = adjacencyList[i].size();
			ig[i + 1] = ig[i] + size;
		}

		jg.resize(ig[2 * slaeSize]);
		int cur = 0;
		for (int i = 1; i < 2 * slaeSize; i++)
			for (int j = 0; j < adjacencyList[i].size(); j++)
				jg[cur++] = adjacencyList[i][j];

		ggl.resize(ig[2 * slaeSize]);
		ggu.resize(ig[2 * slaeSize]);
	}

	//��������� ������� �� ������
	void Matrix::MultiplyAx(const vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;

		for (i = 0; i < n; i++)
		{
			//������ i-�� ������(�������)
			l = ig[i];
			//������ (i+1)-�� ������(�������)
			iend = ig[i + 1];
			//���������� ��������� � i ������(�������)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//�������� �� ���� ��������� i ������ (�������)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggl[l] * a[j];
				result[j] += ggu[l] * a[i];
			}
		}
	}

	//��������� ����������������� ������� �� ������
	void Matrix::MultiplyATx(vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//������ i-�� ������(�������)
			l = ig[i];
			//������ (i+1)-�� ������(�������)
			iend = ig[i + 1];
			//���������� ��������� � i ������(�������)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//�������� �� ���� ��������� i ������ (�������)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggu[l] * a[j];
				result[j] += ggl[l] * a[i];
			}
		}
	}

}

