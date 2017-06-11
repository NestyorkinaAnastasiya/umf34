/*slae.h*/
#pragma once
#include "tests.h"
#include "LU.h"
#include <ctime>
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace tests;
using namespace lu;

namespace slae
{
	class SLAE : private Basis, private GaussIntegration
	{
		// ����������� ������
		int n;
		// ������������ ���������� �������� � ��������
		int maxiter = 10000;
		// �������� ������� ����
		const double eps = 1e-14;
		// �������
		double omega;
		// �����
		Grid grid;
		// ��������� �������� �������
		Tests tests;
		// ���������� �������
		Matrix A;
		// ��������� �������
		// ������� ��������
		array<array<double, CountOfDofs>, CountOfDofs> G;
		// ������� �����
		array<array<double, CountOfDofs>, CountOfDofs> MSigma, MEps;
		// ��������� ������� ������ �����
		array <double, CountOfDofs> locFSin, locFCos;
		// ���������� ������ ������ �����
		vector <double> F;
		// ������ ������������� �������
		vector <double> u;
		// ����� ������� ������ �����
		double normF;

		// ������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		// ������ ��������� ������ ����
		void CalculateM(int elementNumber);
		// ������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		// ������� ���������� �������� � ����������
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		// ������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);

		// ������ ������� ����� ��� ������� �������� 
		// 0 - sin, 1 - cos
		array<array<double, 2>,2> g;
		// ���������� ������ ����� ��� 1��� �������� �������
		void Calculate_g(int formNumber, int orientation, int elNumber);
		// ���������� 1��� �������� ������� ��� ������ ����
		void CalculateBoundaries1ForNode(int node, double weight);
		// ���� ������� �������� �������
		void CalculateBoundaries1(int number);

		// ���������� ������� � �������������
		vector <double> L;
		vector <double> D;
		vector <double> U;

		// ������� �������� �������
		double t;
		// ������ �������
		vector <double> r;
		// ������ ������
		vector <double> z;

		// ���������� ����� �������
		double Norm(const vector<double>& x);
		// ��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		// ��������� ���� �� i-�� �������� �� �������
		void GenerateSLAE();
		// LU-������������
		void LU();
		// ��������������� ������� ��� ��������
		void LYF(const vector<double> C, vector<double>& yl);
		void UXY(const vector<double> C, vector<double>& yu);
		void UTYF(const vector<double> f, vector<double>& y);
		void LTXY(const vector<double> y, vector<double>& x);
		double Rel_Discrepancy();

		lu::LU solverLU;

		// �������� ��� � LU-�������������
		void LULOS();
		void SolveLU();
		void BCG_STAB();
		void LUBCG();
	public:
		SLAE();
		void Solve();
		~SLAE() {};
	};
}
