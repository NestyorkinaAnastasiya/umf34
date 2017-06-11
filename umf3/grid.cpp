/*grid.cpp*/
#include "grid.h"

namespace grid
{
	void Area::Input(FILE *fo)
	{
		fscanf_s(fo, "%d", &leftX);
		fscanf_s(fo, "%d", &rightX);
		fscanf_s(fo, "%d", &lowY);
		fscanf_s(fo, "%d", &upY);

		fscanf_s(fo, "%d", &ku[0]);
		fscanf_s(fo, "%d", &ku[1]);
		fscanf_s(fo, "%d", &ku[2]);
		fscanf_s(fo, "%d", &ku[3]);

		fscanf_s(fo, "%d", &kuForm[0]);
		fscanf_s(fo, "%d", &kuForm[1]);
		fscanf_s(fo, "%d", &kuForm[2]);
		fscanf_s(fo, "%d", &kuForm[3]);
	}

	AreasLines::AreasLines()
	{
		int n;
		double tmp;
		FILE *fo;
		fopen_s(&fo,"AreasLines.txt", "r");
		fscanf_s(fo, "%d", &n);
		x.reserve(n);
		for (int i = 0; i < n; i++)
		{
			fscanf_s(fo, "%lf", &tmp);
			x.push_back(tmp);
		}

		fscanf_s(fo, "%d", &n);
		y.reserve(n);
		for (int j = 0; j < n; j++)
		{
			fscanf_s(fo, "%lf", &tmp);
			y.push_back(tmp);
		}
		fclose(fo);
	}

	//Поиск глобального номера б.ф. dofsNumber в элементе
	bool Element::SearchDof(int dofsNumber)
	{
		for (int i = 0; i < CountOfDofs; i++)
			if (dofsNumber == dof[i]) return true;

		return false;
	}

	Grid::Grid()
	{
		int n;
		FILE *fo;
		fopen_s(&fo,"Areas.txt", "r");
		fscanf_s(fo, "%d", &n);
		areas.resize(n);
		for (int i = 0; i < n; i++)
			areas[i].Input(fo);
		fclose(fo);
	}

	Grid::~Grid() {}

	//Генерация координаты с учётом разбиения на всех интервалах
	void Grid::PartitionСoordinate(vector <double> &x, vector <double> areasLines,
		vector <double> coefficient, vector <int> nIntervals)
	{
		int count;
		//длина интервала
		double l;
		//шаг
		double h;
		//число интервалов
		int nLines = areasLines.size();

		count = 0;
		for (int i = 0; i < nLines - 1; i++)
		{
			x.push_back(areasLines[i]);
			count++;

			//длина интервала
			l = abs(areasLines[i + 1] - areasLines[i]);

			//рассчитываем первый шаг
			//равномерная
			if (abs(1.0 - coefficient[i]) < 1E-14)
				h = l / nIntervals[i];
			else //сгущаем вправо
				if (coefficient[i] < 1)
				{
					h = l * (1.0 - coefficient[i]);
					h /= 1.0 - pow(coefficient[i], nIntervals[i]);
				}
				else //сгущаем влево
				{
					h = l * (coefficient[i] - 1.0);
					h /= coefficient[i] * (pow(coefficient[i], nIntervals[i] - 1.0) - 1.0);
					h += coefficient[i] - 1.0;
				}

			//получаем сетку внутри интервала
			for (int j = 1; j < nIntervals[i]; j++)
			{
				if (j != 1) h *= coefficient[i];
				x.push_back(x[count - 1] + h);
				count++;
			}
		}
		x.push_back(areasLines[nLines - 1]);
	}

	//Добавление узла
	void Grid::PushNode(double x, double y)
	{
		Point p;
		p.y = y;
		p.x = x;
		//кладём узел в конец массива узлов
		nodes.push_back(p);
	}

	//Построение сетки
	void Grid::BuildGrid()
	{
		//для iого интервала и каждой координаты
		//коэффициент разрядки
		vector <double> xCoefficient;
		vector <double> yCoefficient;
		//число подинтервалов
		vector <int> xIntervals;
		vector <int> yIntervals;

		//количество координатных линий по x и y
		int xLines = areasLines.x.size();
		int yLines = areasLines.y.size();

		//геометрические линии разбиения по х и y
		vector <double> xi;
		vector <double> yj;

		//временная переменная для считывания
		int tmp1;
		//временная переменная для считывания
		double tmp2;

		//число узлов и число элементов
		int nNodes, nElements;

		xIntervals.reserve(xLines - 1); xCoefficient.reserve(xLines - 1);
		yIntervals.reserve(yLines - 1); yCoefficient.reserve(yLines - 1);

		FILE *fo;
		fopen_s(&fo,"Intervals.txt", "r");
		//ввод количества интервалов по х и у и коэффициентов разрядки
		for (int i = 0; i < xLines - 1; i++)
		{
			fscanf_s(fo, "%d", &tmp1);
			xIntervals.push_back(tmp1);
			fscanf_s(fo, "%lf", &tmp2);
			xCoefficient.push_back(tmp2);
		}

		for (int j = 0; j < yLines - 1; j++)
		{
			fscanf_s(fo, "%d", &tmp1);
			yIntervals.push_back(tmp1);
			fscanf_s(fo, "%lf", &tmp2);
			yCoefficient.push_back(tmp2);
		}
		fclose(fo);
		nx = 0;
		for (int i = 0; i < xLines - 1; i++)
			//общее количество интервалов по х
			nx += xIntervals[i];
		nx++;

		ny = 0;
		for (int j = 0; j < yLines - 1; j++)
			//общее количество интервалов по у
			ny += yIntervals[j];
		ny++;

		xi.reserve(nx); yj.reserve(ny);

		//построение сеток по х и у
		PartitionСoordinate(xi, areasLines.x, xCoefficient, xIntervals);
		PartitionСoordinate(yj, areasLines.y, yCoefficient, yIntervals);

		xIntervals.clear(); xCoefficient.clear();
		yIntervals.clear(); yCoefficient.clear();


		nElements = (nx - 1) * (ny - 1);
		nNodes = xi.size() * yj.size();

		elements.reserve(nElements);
		nodes.reserve(nNodes);

		//заполняем список узлов
		for (int j = 0; j < yj.size(); j++)
			for (int i = 0; i < xi.size(); i++)
				PushNode(xi[i], yj[j]);
	}

	//Получение глобального номера
	int Grid::GetGlobalNumber(int elementNumber, int localNumber)
	{
		//число интервалов по x
		int nxint = nx - 1;
		//номер горизонтальной линии, которая является нижним ребром элемента
		int nGorLine = elementNumber / nxint;
		//начальный номер узна на ребре
		int nodeOnLine = nGorLine * (nxint + 1);
		//отступ от начального узла, чтобы получить текущий номер
		//левого нижнего узла элемента;
		int offsetX = elementNumber % nxint;
		//номер нижнего левого узла
		int nodeElem = nodeOnLine + offsetX;


		if (localNumber == 0)
			return nodeElem;

		if (localNumber == 1)
			return nodeElem + 1;

		if (localNumber == 2)
			return nodeElem + (nxint + 1);

		if (localNumber == 3)
			return nodeElem + (nxint + 1) + 1;

		return -1;

	}

	//Нахождение номера области, в которой лежит заданная точка
	int Grid::FindArea(double x, double y)
	{
		bool xInterval, yInterval;
		int size = areas.size();
		for (int i = 0; i < size; i++)
		{
			xInterval = x > areasLines.x[areas[i].leftX] && x < areasLines.x[areas[i].rightX];
			yInterval = y > areasLines.y[areas[i].lowY] && y < areasLines.y[areas[i].upY];

			//если точка попадает в iую область, возвращаем номер этой области
			if (xInterval && yInterval) return i;
		}
	}

	//Нахождение соседних конечных элементов
	void Grid::FindNeighbors(int elementNumber)
	{
		int count = 0;
		//количество конечных элементов
		int nElements = elements.size();
		Element element;
		element = elements[elementNumber];

		bool tmp[4] = { false, false, false, false };

		for (int i = 0; count < 4 && i < nElements; i++)
		{
			if (i != elementNumber)
			{
				//сосед по левому ребру
				if (nodes[element.nodes[0]] == nodes[elements[i].nodes[1]] && nodes[element.nodes[2]] == nodes[elements[i].nodes[3]])
				{
					element.neighbors[0] = i;
					tmp[0] = true;
					count++;
				}
				//сосед по правому ребру
				if (nodes[element.nodes[1]] == nodes[elements[i].nodes[0]] && nodes[element.nodes[3]] == nodes[elements[i].nodes[2]])
				{
					element.neighbors[1] = i;
					tmp[1] = true;
					count++;
				}
				//сосед по нижнему ребру
				if (nodes[element.nodes[0]] == nodes[elements[i].nodes[2]] && nodes[element.nodes[1]] == nodes[elements[i].nodes[3]])
				{
					element.neighbors[2] = i;
					tmp[2] = true;
					count++;
				}
				//сосед по верхнему ребру
				if (nodes[element.nodes[2]] == nodes[elements[i].nodes[0]] && nodes[element.nodes[3]] == nodes[elements[i].nodes[1]])
				{
					element.neighbors[3] = i;
					tmp[3] = true;
					count++;
				}
			}
		}

		for (int i = 0; i < 4; i++)
			//отсутствие соседа на соответствующей стороне
			if (tmp[i] == false) element.neighbors[i] = -1;

		elements[elementNumber] = element;
	}

	//Вычисление конечных элементов
	void Grid::ComputeElements()
	{
		Element tmpEl;
		double x, y;

		//вычисляем глобальные номера узлов и базисных функций
		//каждого кэ и заполняем список кэ
		for (int i = 0; i < elements.capacity(); i++)
		{
			for (int j = 0; j < 4; j++)
				tmpEl.nodes[j] = GetGlobalNumber(i, j);

			for (int j = 0; j < CountOfDofs; j++)
				tmpEl.dof[j] = GetGlobalNumber(i, j);

			elements.push_back(tmpEl);
		}

		//находим номер подобласти для каждого кэ
		for (int i = 0; i < elements.size(); i++)
		{

			//берём центральный узел
			x = (nodes[elements[i].nodes[0]].x + nodes[elements[i].nodes[1]].x) / 2;
			y = (nodes[elements[i].nodes[0]].y + nodes[elements[i].nodes[2]].y) / 2;

			elements[i].numberOfArea = FindArea(x, y);
		}

		//находим соседние кэ
		for (int i = 0; i < elements.size(); i++)
			FindNeighbors(i);
	}

	//Формирование массивов краевых условий
	void Grid::FormKU()
	{
		int size = elements.size();
		BoundaryCondition tmp;

		//для каждого элемента находим,какие ку на его границах
		for (int i = 0; i < size; i++)
		{
			bool left[3] = { false, false, false }, right[3] = { false, false, false },
				low[3] = { false, false, false }, up[3] = { false, false, false };
			if (elements[i].neighbors[0] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[0] == 1) left[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[0] == 2) left[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[0] == 3) left[2] = true;
			}
			if (elements[i].neighbors[1] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[1] == 1) right[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[1] == 2) right[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[1] == 3) right[2] = true;
			}
			if (elements[i].neighbors[2] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[2] == 1) low[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[2] == 2) low[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[2] == 3) low[2] = true;
			}
			if (elements[i].neighbors[3] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[3] == 1) up[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[3] == 2) up[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[3] == 3) up[2] = true;
			}

			for (int j = 0; j < 3; j++)
				if (left[j] || right[j] || low[j] || up[j])
				{
					tmp.elem = i;
					if (left[j])
					{
						tmp.edges[0] = 1;
						tmp.formNumber[0] = areas[elements[i].numberOfArea].kuForm[0];
					}
					else
					{
						tmp.edges[0] = 0;
						tmp.formNumber[0] = -1;
					}
					if (right[j])
					{
						tmp.edges[1] = 1;
						tmp.formNumber[1] = areas[elements[i].numberOfArea].kuForm[1];
					}
					else
					{
						tmp.edges[1] = 0;
						tmp.formNumber[1] = -1;
					}
					if (low[j])
					{
						tmp.edges[2] = 1;
						tmp.formNumber[2] = areas[elements[i].numberOfArea].kuForm[2];
					}
					else
					{
						tmp.edges[2] = 0;
						tmp.formNumber[2] = -1;
					}
					if (up[j])
					{
						tmp.edges[3] = 1;
						tmp.formNumber[3] = areas[elements[i].numberOfArea].kuForm[3];
					}
					else
					{
						tmp.edges[3] = 0;
						tmp.formNumber[3] = -1;
					}

					ku[j].push_back(tmp);
				}
		}
	}

	//Образование множества кэ
	void Grid::DoPartition()
	{
		//Построение сетки
		BuildGrid();
		//Вычисление конечных элементов
		ComputeElements();
		//Формирование массивов краевых условий
		FormKU();
	}

	//Формирование списка элементов, содержащих глобальный номер б.ф.
	//равный dofsNumber
	void Grid::SearchElements(int dofsNumber, vector <int> &elList)
	{
		int count;
		int size = elements.size();
		elList.reserve(4);

		count = 0;
		for (int i = 0; i < size && count < 4; i++)
		{
			if (elements[i].SearchDof(dofsNumber))
			{
				count++;
				elList.push_back(i);
			}
		}
	}
}
