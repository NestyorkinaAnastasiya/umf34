/*slae.cpp*/
#include "slae.h"

namespace slae
{
	SLAE::SLAE()
	{
		//Построение сетки
		grid.DoPartition();
		//Размерность задачи соответствует общему числу базисных функций
		n = grid.nx * grid.ny;

		//Генерация портрета матрицы и её инициализация
		A.CreatePortret(n, grid);
		n *= 2;
		F.resize(n);
		u.resize(n);
		r.resize(n);
		z.resize(n);
	}

	void SLAE::Solve()
	{
		GenerateSLAE();
		vector<double> x(n), ut(n);
		LULOS();
		SolveLU();
		BCG_STAB();
	}

	//Сборка локальных матриц жёсткости
	void SLAE::CalculateG(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double g1, g2, ksi, etta, x_, y_, lambda,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			hx2 = hx * hx, hy2 = hy * hy,
			jacobian = hx * hy / 4.0;

		for (int i = 0; i < CountOfDofs; i++)
		{	
			for (int j = i; j < CountOfDofs; j++)
			{	
				g1 = 0; g2 = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;
					lambda = tests.Lambda(x_, y_);

					g1 += gaussWeights[k] * dphiksi[i](ksi, etta) * dphiksi[j](ksi, etta) * lambda;
					g2 += gaussWeights[k] * dphietta[i](ksi, etta) * dphietta[j](ksi, etta) * lambda;
				}
				G[i][j] = g1 * jacobian / hx2 + g2 * jacobian / hy2;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < CountOfDofs; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];
	}

	//Сборка локальных матриц масс
	void SLAE::CalculateM(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double  g1, g2, ksi, etta, sigma, eps, x_, y_,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;

		omega = tests.Omega();
		
		for (int i = 0; i < CountOfDofs; i++)
		{	
			for (int j = i; j < CountOfDofs; j++)
			{	
				g1 = g2 = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;


					g1 += tests.Sigma(x_, y_) *  gaussWeights[k] * phi[i](ksi, etta) * phi[j](ksi, etta);
					g2 += tests.Eps(x_, y_) *  gaussWeights[k] * phi[i](ksi, etta) * phi[j](ksi, etta);
				}
				MSigma[i][j] = g1 * jacobian * omega;
				MEps[i][j] = g2 * jacobian * omega * omega;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < CountOfDofs; i++)
			for (int j = 0; j < i; j++)
			{
				MSigma[i][j] = MSigma[j][i];
				MEps[i][j] = MEps[j][i];
			}
				
	}

	//Сборка локальных правых частей
	void SLAE::CalculateLocalF(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double fi, ksi, etta, x_, y_, ui,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;
		//интегрирование(Гаусс 5)
		for (int i = 0; i < CountOfDofs; i++)
		{	
			locFSin[i] = 0;
			locFCos[i] = 0;
			
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;

					locFSin[i] += tests.FiSin(x_, y_) * gaussWeights[k] * phi[i](ksi, etta);
					locFCos[i] += tests.FiCos(x_, y_) * gaussWeights[k] * phi[i](ksi, etta);
				}
			locFSin[i] *= jacobian ;
			locFCos[i] *= jacobian;
		}
	}

	//Добавка локального элемента в глобальный
	void SLAE::AddElementToGlobalMatrix(Matrix &B, int i, int j, double element)
	{
		int id;
		bool flag;

		if (i == j)
			B.di[i] += element;
		else
		{
			if (i < j)
			{
				flag = false;
				for (id = B.ig[j]; !flag && id < B.ig[j + 1]; id++)
					if (B.jg[id] == i) flag = true;

				if (flag) B.ggu[id - 1] += element;
			}
			else
			{
				flag = false;
				for (id = B.ig[i]; !flag && id < B.ig[i + 1]; id++)
					if (B.jg[id] == j) flag = true;

				if (flag) B.ggl[id - 1] += element;
			}
		}
	}

	// Сборка локальных матриц(векторов) и добавление в глобальные
	void SLAE::CalculateLocals(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		int ki, kj;

		// Вычисление локальных матриц
		CalculateG(elementNumber);
		CalculateM(elementNumber);
		CalculateLocalF(elementNumber);

		for (int i = 0; i < CountOfDofs; i++)
		{
			ki = element.dof[i];
			for (int j = 0; j < CountOfDofs; j++)
			{
				kj = element.dof[j];
				// Добавка в глобальную матрицу А
				// G*qsin
				AddElementToGlobalMatrix(A, 2 * ki, 2 * kj, G[i][j] - MEps[i][j]);
				// -M*qcos
				AddElementToGlobalMatrix(A, 2 * ki, 2 * kj + 1, -MSigma[i][j]);
				// M*qsin
				AddElementToGlobalMatrix(A, 2 * ki + 1, 2 * kj, MSigma[i][j]);
				// G*qcos
				AddElementToGlobalMatrix(A, 2 * ki + 1, 2 * kj + 1, G[i][j] - MEps[i][j]);
			}
			//добавка в глобальную правую часть
			F[2 * ki] += locFSin[i];
			F[2 * ki + 1] += locFCos[i];
		}
	}

	//Генерация СЛАУ
	void SLAE::GenerateSLAE()
	{
		for (int i = 0; i < n; i++)
			A.di[i] = F[i] = 0;
		for (int i = 0; i < A.ggl.size(); i++)
			A.ggl[i] = A.ggu[i] = 0;

		//Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		//Учёт краевых условий
		for (int i = 0; i < grid.ku[0].size(); i++)
			CalculateBoundaries1(i);

		normF = Norm(F);
	}

	//Нахождение правой части для 1ого краевого условия
	void SLAE::Calculate_g(int formNumber, int orientation, int elNumber)
	{
		Element element = grid.elements[elNumber];

		switch (orientation)
		{
		//левое ребро
		case 0:
		{
			double x = grid.nodes[element.nodes[0]].x;
			double y1 = grid.nodes[element.nodes[0]].y;
			double y2= grid.nodes[element.nodes[2]].y;

			g[0][0] = tests.UgSin(x, y1);
			g[0][1] = tests.UgSin(x, y2);

			g[1][0] = tests.UgCos(x, y1);
			g[1][1] = tests.UgCos(x, y2);
		}
		break;
		//правое ребро
		case 1:
		{
			double x = grid.nodes[element.nodes[1]].x;
			double y1 = grid.nodes[element.nodes[1]].y;
			double y2 = grid.nodes[element.nodes[3]].y;

			g[0][0] = tests.UgSin(x, y1);
			g[0][1] = tests.UgSin(x, y2);

			g[1][0] = tests.UgCos(x, y1);
			g[1][1] = tests.UgCos(x, y2);
		}
		break;
		//нижнее ребро
		case 2:
		{
			double y = grid.nodes[element.nodes[0]].y;
			double x1 = grid.nodes[element.nodes[0]].x;
			double x2 = grid.nodes[element.nodes[1]].x;
			g[0][0] = tests.UgSin(x1, y);
			g[0][1] = tests.UgSin(x2, y);

			g[1][0] = tests.UgCos(x1, y);
			g[1][1] = tests.UgCos(x2, y);
		}
		break;
		//верхнее ребро
		case 3:
		{
			double y = grid.nodes[element.nodes[2]].y;
			double x1 = grid.nodes[element.nodes[2]].x;
			double x2 = grid.nodes[element.nodes[3]].x;
			g[0][0] = tests.UgSin(x1, y);
			g[0][1] = tests.UgSin(x2, y);

			g[1][0] = tests.UgCos(x1, y);
			g[1][1] = tests.UgCos(x2, y);
		}
		break;
		default:; break;
		}
	}

	//Вычисление 1ого краевого условия для одного узла
	void SLAE::CalculateBoundaries1ForNode(int node, double weight)
	{
		int id;
		F[node] = weight;
		A.di[node] = 1;

		for (int j = 0; j < n; j++)
			if (node < j)
			{
				bool flag = false;
				for (id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
					if (A.jg[id] ==  node) flag = true;
				if (flag) A.ggu[id - 1] = 0.0;
			}
			else
			{
				bool flag = false;
				for (id = A.ig[node]; !flag && id <= A.ig[node + 1] - 1; id++)
					if (A.jg[id] == j) flag = true;
				if (flag) A.ggl[id - 1] = 0.0;
			}
	}

	//Учёт первого краевого условия
	void SLAE::CalculateBoundaries1(int number)
	{
		Element element = grid.elements[grid.ku[0][number].elem];
		// Левая граница
		if (grid.ku[0][number].edges[0] == 1)
		{
			int indexes[2] = { element.dof[0], element.dof[2]};
			Calculate_g(grid.ku[0][number].formNumber[0], 0, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
			{ 
				CalculateBoundaries1ForNode(2 * indexes[i], g[0][i]);
				CalculateBoundaries1ForNode(2 * indexes[i] + 1, g[1][i]);
			}

		}
		// Правая граница
		if (grid.ku[0][number].edges[1] == 1)
		{
			int indexes[2] = { element.dof[1], element.dof[3]};
			Calculate_g(grid.ku[0][number].formNumber[1], 1, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
			{
				CalculateBoundaries1ForNode(2 * indexes[i], g[0][i]);
				CalculateBoundaries1ForNode(2 * indexes[i] + 1, g[1][i]);
			}
		}
		// Нижняя граница
		if (grid.ku[0][number].edges[2] == 1)
		{
			int indexes[2] = { element.dof[0], element.dof[1]};
			Calculate_g(grid.ku[0][number].formNumber[2], 2, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
			{
				CalculateBoundaries1ForNode(2 * indexes[i], g[0][i]);
				CalculateBoundaries1ForNode(2 * indexes[i] + 1, g[1][i]);
			}
		}
		// Верхняя граница
		if (grid.ku[0][number].edges[3] == 1)
		{
			int indexes[2] = { element.dof[2], element.dof[3] };
			Calculate_g(grid.ku[0][number].formNumber[3], 3, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
			{
				CalculateBoundaries1ForNode(2 * indexes[i], g[0][i]);
				CalculateBoundaries1ForNode(2 * indexes[i] + 1, g[1][i]);
			}
		}
	}

	//Вычисление нормы вектора
	double SLAE::Norm(const vector<double> &x)
	{
		double norm = 0;
		int size = x.size();

		for (int i = 0; i < size; i++)
			norm += x[i] * x[i];

		return sqrt(norm);
	}

	//Скалярное произведение векторов
	double SLAE::Scalar(const vector<double> &x, const vector<double> &y)
	{
		double sum = 0;
		int size = x.size();
		for (int i = 0; i < size; i++)
			sum += x[i] * y[i];

		return sum;
	}

	//LU-факторизация
	void SLAE::LU()
	{
		int i, i0, j0, iend, num, ki, kj, jend;
		double suml, sumu, sumdg;

		L.resize(A.ggl.size());
		L = A.ggl;
		U.resize(A.ggu.size());
		U = A.ggu;
		D.resize(A.di.size());
		D = A.di;

		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			for (num = i0, sumdg = 0; num < iend; num++)
			{
				j0 = A.ig[A.jg[num]];
				jend = A.ig[A.jg[num] + 1];
				ki = i0;
				kj = j0;
				//для num учитываются все предыдущие элементы
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++)
				{
					for (int m = kj; m < jend; m++)
						//ищем соответствующие ненулевые элементы для умножения
						if (A.jg[ki] == A.jg[m])
						{
							suml += L[ki] * U[m];
							sumu += L[m] * U[ki];
						}
				}
				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[A.jg[num]];
				//умножаются симметричные элементы
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}
	
	void SLAE::LYF(const vector <double> f, vector <double> &y)
	{
		//i0-адрес начала строки, iend-адрес конца строки
		int i, i0, iend;
		double sum;
		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i]; iend = A.ig[i + 1];

			y[i] = f[i];

			for (i0, sum = 0; i0 < iend; i0++)
				y[i] -= y[A.jg[i0]] * L[i0];
			y[i] /= D[i];
		}
	}

	void SLAE::UXY(const vector <double> y, vector <double> &x)
	{
		int i, i0, iend;

		for (i = 0; i < n; i++)
			x[i] = 0.0;
		//проход по столбцам с конца
		for (i = n - 1; i >= 0; i--)
		{
			x[i] += y[i];

			i0 = A.ig[i]; iend = A.ig[i + 1]; iend--;
			//идём по столбцу с конца
			for (; iend >= i0; iend--)
				x[A.jg[iend]] -= x[i] * U[iend];
		}
	}

	double SLAE::Rel_Discrepancy()
	{
		double dis1, dis2;
		dis1 = Scalar(r, r);
		dis2 = Scalar(F, F);
		double dis = dis1 / dis2;
		return sqrt(dis);
	}

	void SLAE::LULOS()
	{
		double a, b, pp, dis, rr;
		int i,k;
		vector <double> Ax(n), C(n), p(n);
		FILE *fo;
		fopen_s(&fo, "LULOS_result.txt", "w");
		for (auto &i : u) i = 0;
		unsigned int beginTime = clock();
		LU();
		//Ax0
		A.MultiplyAx(u, Ax);
		//f-Ax0
		for (i = 0; i < n; i++)
			r[i] = F[i] - Ax[i];

		//r0=L^(-1)(f-Ax0)
		LYF(r, r);

		//z0=U^(-1)r0->r0=Uz0
		UXY(r, z);

		//p0=L^(-1)Az0
		A.MultiplyAx(z, Ax);//Az0
		LYF(Ax, p);

		rr = Scalar(r, r);
		dis = Scalar(r, r) / rr;
		dis = sqrt(dis);
		k = 0;

		for (; dis > eps && k <= maxiter; k++)
		{
			//Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			//Xk, Rk
			for (i = 0; i < n; i++)
			{
				u[i] = u[i] + a*z[i];
				r[i] = r[i] - a*p[i];
			}

			//UY=rk->Y=U^(-1)rk
			UXY(r, C);
			//AU^(-1)rk=Ax
			A.MultiplyAx(C, Ax);
			//L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
			LYF(Ax, Ax);

			//bk
			b = -Scalar(p, Ax) / pp;

			//zk=U^(-1)rk+bkz[k-1]
			//pk
			for (i = 0; i < n; i++)
			{
				z[i] = C[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}
			dis = Scalar(r, r) / rr;
			dis = sqrt(dis);
		}

		unsigned int endTime = clock();
		dis = Rel_Discrepancy();

		fprintf(fo, "time: %d мс\n", endTime - beginTime);
		fprintf(fo, "iter: %d\n", k);
		fprintf(fo, "res: %e\n", dis);
		
		for (int i = 0; i < n; i++)
		{
			fprintf(fo, "%.14f\t", u[i]);
			fprintf(fo, "%.14f\n", u[i + 1]);
			i++;
		}
		fclose(fo);
	}

	void SLAE::SolveLU()
	{
		FILE *fo;
		fopen_s(&fo, "LU_result.txt", "w");
		for (auto &i : u) i = 0;
		unsigned int beginTime = clock();
		solverLU.Convert(A.ig, A.jg, A.di, F, A.ggl, A.ggu, n);
		solverLU.Solve(u);
		unsigned int endTime = clock();
		fprintf(fo, "time: %d мс\n", endTime - beginTime);
		
		for (int i = 0; i < n; i++)
		{
			fprintf(fo, "%.14f\t", u[i]);
			fprintf(fo, "%.14f\n", u[i+1]);
			i++;
		}
		fclose(fo);
	}

	void SLAE::BCG_STAB()
	{
		double alpha, betta, gamma, rr0, dis, pr;
		int i, k;
		vector <double> r0(n), rPrev(n), tmpZ(n), tmpP(n), p(n);
		FILE *fo;
		fopen_s(&fo, "BCG_STAB_result.txt", "w");
		for (auto &i : u) i = 0;

		unsigned int beginTime = clock();
		LU();
		
		//Ax0
		A.MultiplyAx(u, r0);
		
		//f-Ax0
		for (i = 0; i < n; i++)
			r0[i] = F[i] - r0[i];

		LYF(r0, r0);
		UXY(r0, z);
		k = 0;
		dis = Norm(r0);
		r = r0;
		for (; dis > eps && k <= maxiter; k++)
		{
			rPrev = r;
			//UY=zk-1->Y=U^(-1)zk-1
			UXY(z, tmpZ);
			//AU^(-1)zk-1
			A.MultiplyAx(tmpZ, tmpZ);
			//L^(-1)AU^(-1)zk-1=Y2->L^(-1)B=Y2->LY2=B
			LYF(tmpZ, tmpZ);

			rr0 = Scalar(rPrev, r0);

			// ak = (r,r0)/ ( L^(-1)AU^(-1)zk,r0)
			alpha = rr0 / Scalar(r0, tmpZ);

			// pk = r - ak*L^(-1)AU^(-1)zk
			for (i = 0; i < n; i++)
				p[i] = rPrev[i] - alpha*tmpZ[i];

			UXY(p, tmpP);
			A.MultiplyAx(tmpP, tmpP);
			LYF(tmpP, tmpP);

			// gk = (L^(-1)AU^(-1)pk,pk)/ ( L^(-1)AU^(-1)pk,L^(-1)AU^(-1)pk)
			gamma = Scalar(p, tmpP) / Scalar(tmpP, tmpP);
						
			for (i = 0; i < n; i++)
			{
				//u =u + ak * zk + gk * pk;
				u[i] = u[i] + alpha*z[i] + gamma*p[i];
				//r = pk - gk * L^(-1)AU^(-1)pk;
				r[i] = p[i] - gamma*tmpP[i];
			}

			betta = alpha * Scalar(r, r0) / gamma / rr0;

			for (i = 0; i < n; i++)
				z[i] = r[i] + betta*z[i] - betta*gamma*tmpZ[i];

			dis = Norm(r);
		}

		UXY(u, u);
		unsigned int endTime = clock();

		fprintf(fo, "time: %d мс\n", endTime - beginTime);
		fprintf(fo, "iter: %d\n", k);
		fprintf(fo, "res: %e\n", dis);
		for (int i = 0; i < n; i++)
		{
			fprintf(fo, "%.14f\t", u[i]);
			fprintf(fo, "%.14f\n", u[i + 1]);
			i++;
		}
		fclose(fo);
	}

}