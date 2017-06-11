/*basis.h*/
#pragma once
#include "grid.h"

using namespace grid;
namespace basis
{
	const int nFunc = 4;
	const int nFunc1D = 2;

	struct Basis
	{
		//Указатели на функции вычисления базисных функций в точке
		array<function<double(double, double)>, nFunc> phi;
		//Указатели на функции вычисления d/dksi базисных функций в точке
		array<function<double(double, double)>, nFunc> dphiksi;
		//Указатели на функции вычисления d/detta базисных функций в точке
		array<function<double(double, double)>, nFunc> dphietta;
		Basis();
	};
}