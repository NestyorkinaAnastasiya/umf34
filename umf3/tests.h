/*tests.h*/
#include "integration.h"

namespace tests
{
	struct Tests
	{
		int test;
		double Lambda(double x, double y);
		double Sigma(double x, double y);
		double Eps(double x, double y);
		double Omega();
		double UgSin(double x, double y);
		double UgCos(double x, double y);
		double FiSin(double x, double y);
		double FiCos(double x, double y);
		Tests();
	};
}
