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
		//Формируем новый массив ig
		ig[0] = 0;
		for (int i = 1; i <= n; i++)
		{
			//Количество не нулевых элементов в строке
			int k = _ig[i] - _ig[i - 1];
			//Если такие есть
			if (k > 0) 
			{	// Вычисляем сколько всего элементов в строке, как в профильной
				int total_n = i - 1 - _jg[_ig[i-1]];
				ig[i] = ig[i - 1] + total_n;
			}
			else
				ig[i] = ig[i - 1];
		}
		U.resize(ig[n]);
		L.resize(ig[n]);

		// Формируем новые gl и gu

		for (int i = 0; i < n; i++) 
		{	// Начало строки
			int start = ig[i];
			// Конец
			int end = ig[i + 1];
			// Номер текущей колонки
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
		//построение LU(sq) разложения
		Decomposition();

		//Решение системы Ly=b (прямой проход) 
		LYF();

		//Решение системы Ux=y (обратный проход)
		UXY();
		u = x;
	}

	void LU::Decomposition()
	{	
		int i, j, i0, j0, iend, jend, kol_i, kol_j, num, first_i, first_j, ki, kj, dif;
		//суммы для элементов в L,U
		double suml, sumu, sumdg;

		for (i = 0; i < n; i++)
		{	//начало i строки (столбца)
			i0 = ig[i];
			//начало i+1 строки (столбца)
			iend = ig[i + 1];
			//количество элементов в i строке (столбце)
			kol_i = iend - i0;
			//номер первого ненулевого элемента i строки (столбца)
			first_i = i - kol_i;
			j = first_i;
			//идём по всем ненулевым элементам строки (столбца)
			for (num = i0, sumdg = 0; j<i; num++, j++)
			{
				j0 = ig[j];
				jend = ig[j + 1];
				kol_j = jend - j0;
				first_j = j - kol_j;
				ki = i0;
				kj = j0;
				//разница в профилях
				//смещение по профилю для корректного умножения
				dif = first_j - first_i;
				//ki - номер первого умножаемого элемента в L,если идём по строке
				if (dif > 0) ki += dif;
				//kj - номер первого умножаемого элемента в U,если идём по строке
				else kj -= dif;
				//для num учитываются все предыдущие элементы
				for (suml = 0, sumu = 0; ki<num; ki++, kj++) 
				{
					suml += L[ki] * U[kj];
					//считается сумма для симметричного элемента из U
					sumu += L[kj] * U[ki];
				}

				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[j];
				//умножаются симметричные элементы
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}

	void LU::LYF()
	{	//// i0 - адрес начала строки
		// iend - адрес конца строки
		// ikol - количество элементов в строке
		// beg - номер столбца первого ненулевого элемента строки
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
		//проход по столбцам с конца
		for (int i = n - 1; i >= 0; i--)
		{
			x[i] += y[i];
			i0 = ig[i]; iend = ig[i + 1];
			ikol = iend - i0;
			beg = i - ikol;
			//идём по столбцу с конца
			for (int k = iend - 1, j = i - 1; j >= beg; j--, k--)
				x[j] -= x[i] * U[k];

		}

	}
}