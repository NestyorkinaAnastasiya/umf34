/*tests.cpp*/
#include "tests.h"

namespace tests
{
	double Tests::Lambda(double x, double y)
	{
		switch (test)
		{
		case 1: return 1;
		case 2: return 1;

		case 3: return 100;
		case 4: return 100;
		case 5: return 100;
		case 6: return 100;
		case 7: return 100;

		case 8: return 100;
		case 9: return 1e+05;

		case 10: return 100;
		case 11: return 100;
		case 12: return 100;

		case 13: return 100;
		case 14: return 100;
		case 15: return 100;

		default: return 0;
		}

	}
	double Tests::Sigma(double x, double y)
	{
		switch (test)
		{
		case 1: return 1;
		case 2: return 1;

		case 3: return 1e+04;
		case 4: return 1e+04;
		case 5: return 1e+04;
		case 6: return 1e+04;
		case 7: return 1e+04;

		case 8: return 1e+04;
		case 9: return 1e+04;

		case 10: return 0;
		case 11: return 1e+04;
		case 12: return 1e+08;

		case 13: return 1e+04;
		case 14: return 1e+04;
		case 15: return 1e+04;

		default: return 0;
		}
	}
	double Tests::Eps(double x, double y)
	{
		switch (test)
		{
		case 1: return 1;
		case 2: return 1;

		case 3: return 8.81e-11;
		case 4: return 8.81e-11;
		case 5: return 8.81e-11;
		case 6: return 8.81e-11;
		case 7: return 8.81e-11;

		case 8: return 8.81e-11;
		case 9: return 8.81e-11;

		case 10: return 8.81e-11;
		case 11: return 8.81e-11;
		case 12: return 8.81e-11;

		case 13: return 8.81e-12;
		case 14: return 8.81e-11;
		case 15: return 1e-10;

		default: return 0;
		}
	}
	double Tests::Omega()
	{
		switch (test)
		{
		case 1: return 1;
		case 2: return 1;

		case 3: return 1e-03;
		case 4: return 0;
		case 5: return 1e+03;
		case 6: return 1e+06;
		case 7: return 1e+09;

		case 8: return 1e+03;
		case 9: return 1e+03;

		case 10: return 1e+03;
		case 11: return 1e+03;
		case 12: return 1e+03;

		case 13: return 1e+03;
		case 14: return 1e+03;
		case 15: return 1e+03;

		default: return 0;
		}
	}

	double Tests::UgSin(double x, double y)
	{
		switch (test)
		{
		case 1: return x + y;
		case 2: return pow(x, 2) + pow(y, 2);
		default: return pow(x, 2) + pow(y, 2);
		}

	}
	double Tests::UgCos(double x, double y)
	{
		switch (test)
		{
		case 1: return x - y;
		case 2: return pow(x, 2) - pow(y, 2);
		default: return pow(x, 2) - pow(y, 2);
		}

	}
	double Tests::FiSin(double x, double y)
	{
		switch (test)
		{
		case 1: return -2 * x;
		case 2: return -4 - 2 * pow(x, 2);
		default: return -4 * Lambda(x, y) - pow(Omega(), 2) *Eps(x, y)*(pow(x, 2) + pow(y, 2)) - Omega()*Sigma(x, y)*(pow(x, 2) - pow(y, 2));
		}
	}

	double Tests::FiCos(double x, double y)
	{
		switch (test)
		{
		case 1: return 2 * y;
		case 2: return 2 * pow(y, 2);
		default: return -pow(Omega(), 2) *Eps(x, y)*(pow(x, 2) - pow(y, 2)) + Omega()*Sigma(x, y)*(pow(x, 2) + pow(y, 2));
		}
		return 0;
	}

	Tests::Tests()
	{
		ifstream fo;
		fo.open("number_of_test.txt");
		fo >> test;
	}
}