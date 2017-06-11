/*grid.h*/
#pragma once
#include <stdio.h>
#include <vector>
#include <fstream>
#include <functional>
#include <array>
using namespace std;

namespace grid
{
	const int CountOfDofs = 4;
	//Нумерация идёт следующим образом: 
	//слева,справа, снизу, сверху
	//кэ - конечный элемент
	//б.ф. - базисная функция

	//Область
	struct Area
	{
		//x-координата левой границы
		int leftX;
		//x-координата правой границы
		int rightX;
		//у-координата нижней границы
		int lowY;
		//у-координата верхней границы
		int upY;
		//1 - 1 рода
		//2 - 2 рода
		//3 - 3 рода
		//-1 - нет краевых условий
		array<int, 4> ku;
		//Cпособ расчёта краевых условий
		//-1 - нет краевых условий
		array<int, 4> kuForm;
		void Input(FILE *fo);
	};

	//Координатные линии областей
	struct AreasLines
	{
		//x,y - массивы значений всех вертикальных 
		//и горизонтальных границ подоблостей Wi

		vector <double> x;
		vector <double> y;
		AreasLines();
	};

	//Конечный элемент
	struct Element
	{
		//Узлы
		array<int, 4> nodes;
		//Степени свободы
		array<int, CountOfDofs> dof;
		//Номер области
		int numberOfArea;
		//Соседние элементы
		array<int, 4> neighbors;

		Element& operator=(Element element)
		{
			for (int i = 0; i < 4; i++)
				nodes[i] = element.nodes[i];
			for (int i = 0; i < CountOfDofs; i++)
				dof[i] = element.dof[i];
			numberOfArea = element.numberOfArea;
			for (int i = 0; i < 4; i++)
				neighbors[i] = element.neighbors[i];
			return *this;
		}
		//Поиск глобального номера б.ф. dofsNumber в элементе
		bool SearchDof(int dofsNumber);
	};

	//Краевые условия для конечного элемента
	struct BoundaryCondition
	{
		//Номер элемента
		int elem;
		//Краевое условие(на ребре): 1 - есть, 0 - нет
		array<int, 4> edges;
		//Номер формулы для расчёта
		array<int, 4> formNumber;
	};

	struct Point
	{
		double x;
		double y;
		Point() {};
		~Point() {};
		Point(double xx, double yy)
		{
			x = xx;
			y = yy;
		}

		bool operator==(Point point)
		{
			if (point.x == x && point.y == y)
				return true;
			else
				return false;
		}

	};

	//Разбиение общей области
	class Grid
	{
		//Генерация координаты с учётом разбиения на всех интервалах
		void PartitionСoordinate(vector <double> &ci, vector <double> areasLines,
			vector <double> coefficient, vector <int> nIntervals);
		//Добавление узла
		void PushNode(double x, double y);
		//Построение сетки
		void BuildGrid();

		//Получение глобального номера
		int GetGlobalNumber(int elementNumber, int localNumber);
		
		//Нахождение номера области, в которой лежит заданная точка
		int FindArea(double x, double y);
		//Нахождение соседних конечных элементов
		void FindNeighbors(int elementNumber);
		//Вычисление конечных элементов
		void ComputeElements();

		//Формирование массивов краевых условий
		void FormKU();

	public:
		Grid();
		//Массив областей
		vector <Area> areas;
		//Координатные линии
		AreasLines areasLines;
		//Массив конечных элементов
		vector <Element> elements;
		//Массив узлов
		vector <Point> nodes;
		//Массив краевых условий
		array<vector <BoundaryCondition>, 3> ku;
		//Общее количество интервалов по x и по y
		//после построения сетки
		int nx, ny;
		//Получение глобального номера для базисных функций
		int GetGlobalFuncNumber(int elNum, int localFuncNum);

		//Образование множества кэ
		void DoPartition();
		//Формирование списка элементов, содержащих глобальный номер б.ф.
		//равный dofsNumber
		void SearchElements(int dofsNumber, vector<int>& elList);

		~Grid();
	};

}