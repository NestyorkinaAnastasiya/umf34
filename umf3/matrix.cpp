/*matrix.cpp*/
#include "matrix.h"

namespace matrix
{
	//Генерация портрета матрицы
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
		//// Сортировка списка смежности по возрастанию
		//sort(adjacencyList.begin(), adjacencyList.end());

		/* Сортировка */
		double tmp;
		int elem, size;
		for (int elem = 0; elem < 2 * slaeSize; elem++)
		{
			size = adjacencyList[elem].size();
			// i - номер прохода
			for (int i = 0; i < size - 1; ++i) 
			{	// внутренний цикл прохода
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
		
		/* Инициализация матрицы */
		Initialize(slaeSize);
}

	//Инициализация матрицы после генерации портрета
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

	//Умножение матрицы на вектор
	void Matrix::MultiplyAx(const vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;

		for (i = 0; i < n; i++)
		{
			//начало i-ой строки(столбца)
			l = ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggl[l] * a[j];
				result[j] += ggu[l] * a[i];
			}
		}
	}

	//Умножение транспонированной матрицы на вектор
	void Matrix::MultiplyATx(vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//начало i-ой строки(столбца)
			l = ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggu[l] * a[j];
				result[j] += ggl[l] * a[i];
			}
		}
	}

}

