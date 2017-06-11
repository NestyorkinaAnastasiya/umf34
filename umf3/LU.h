#pragma once
#include <stdio.h>
#include <vector>
using namespace std;
namespace lu
{
	class LU
	{
	public:
		LU();
		~LU();
		int n;
		vector<int> ig;
		vector<double> y, x;
		vector<double>  L, U, D, F;
		//конвертация из строчно-столбцового в профильный
		void Convert(vector<int> _ig, vector<int> _jg, vector<double> _di, vector<double> _b, vector<double> _gl, vector<double> _gu, int size);
		void Solve(vector<double> &u);
	private:
		void Decomposition();
		void LYF();
		void UXY();
	};
}