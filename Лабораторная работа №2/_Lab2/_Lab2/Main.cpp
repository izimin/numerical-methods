#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>
#include <cmath>
#include "linalg.h"

using namespace std;
using namespace alglib;

int const N = 4;
typedef vector<vector<double>> myMatrix;

// Перегрузка оператор >> для вектора
istream &operator >> (istream &in, vector<double> &v)
{
	double num;
	for (int i = 0; i < N; i++) 
	{
		in >> num;
		v.push_back(num);
	}
	return in;
}

// Перегрузка оператор << для вектора
ostream &operator << (ostream &out, vector<double> &v)
{
	for (int i = 0; i < N; i++)
		// 8 знаков после запятой и ширина столбца 15 символов
		out << fixed << setprecision(7) << setw(13) << v[i];
	out << endl;
	return out;
}

// Умножение вектора на число (переопределение оператора *)
vector<double> operator * (vector<double> v, double num)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v[i] * num);
	return res;
}

// Деление вектора на число (переопределение оператора /)
vector<double> operator / (vector<double> v, double num)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v[i] / num);
	return res;
}

// Разность 2-х векторов поэлементно (переопределение оператора -)
vector<double> operator	- (vector<double> v1, vector<double> v2)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v1[i] - v2[i]);
	return res;
}

// Создаёт пустую матрицу
myMatrix EmptyMatrix()
{
	vector<double> v(N);
	myMatrix emp;
	for (int i = 0; i < N; i++)
		emp.push_back(v);
	return emp;
}

// Перегрузка оператора >> для матрицы
istream &operator >> (istream &in, myMatrix &v)
{
	vector<double> str;
	for (int i = 0; i < N; i++)
	{
		in >> str;
		v.push_back(str);
		str.clear();
	}
	return in;
}

// Перегрузка оператора << для матрицы
ostream &operator << (ostream &out, myMatrix &v)
{
	for (int i = 0; i < N; i++)
		// 8 знаков после запятой и ширина столбца 15 символов
		out << fixed << setprecision(7) << setw(13) << v[i];
	return out;
}

// Перегрузка оператора * для перемножения двух матриц
myMatrix operator * (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			for (int k = 0; k < M1[0].size(); k++)
				res[i][j] += M1[i][k] * M2[k][j];
	return res;
}

//  Перегрузка оператора * для умножения матрицы на вектор 
vector<double> operator * (myMatrix &M, vector<double> &v)
{
	vector<double> res(v.size());
	for (int i = 0; i < M.size(); i++)
		for (int j = 0; j < v.size(); j++)
			res[i] += v[j] * M[i][j];			
	return res;
}

// Перегрузка оператора - для перемножения двух матриц
myMatrix operator - (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			res[i][j] += M1[i][j] - M2[i][j
];
	return res;
}

// Функция поиска индекса строки в которой лежит максимальный элемент нужного столбцы
int IndexMaxNumInColumn(myMatrix M, int columnIndex)
{
	int ind = columnIndex;
	for (int i = columnIndex; i < N; i++)
		if (fabs(M[i][columnIndex]) > fabs(M[ind][columnIndex]))
			ind = i;
	if (M[ind][columnIndex] == 0) // Получили 0 на диагонали
		return -1;
	return ind;
}

// Ф-я, которая меняет местами строки 
void Swap(myMatrix &M, int index1, int index2)
{
	auto tmp = M[index1];
	M[index1] = M[index2];
	M[index2] = tmp;
}

// Создаём матрицу U
myMatrix CreateMatrixU(myMatrix &A, myMatrix &P, ofstream &fout, int &rang, double &k)
{
	auto U = A;
	fout << "Построение матрицы U: " << endl;
	int ind;

	// Пробегаем по всем строкам
	for (int i = 0; i < N; i++)
	{
		fout << "\nШАГ " << i + 1 << endl << endl;

		// Ищем индекс максимального элемента в столбце
		ind = IndexMaxNumInColumn(U, i);

		fout << "Максимальный элемент в столбце " << i + 1 << ": " << fabs(U[ind][i]) << endl;

		// Значит весь столбец (начиная со строки i и до N) состоит из 0, а значит ранг равен пройденному числу строк
		if (ind == -1)
		{
			rang = i;
			return U;    // Что-то вернуто-то надо
		}
		
		// Если найденный индекс не явл. индексом тек. строки, то меняем строки
		if (ind != i) 
		{
			k *= -1;
			Swap(U, i, ind);
			Swap(A, i, ind);
			Swap(P, i, ind);
			fout << "Меняем местами " << i + 1 << "-ю и " << ind + 1 << "-ю строки \n" << U;
		}
		else fout << "Перестановка строк не требуется \n";

		fout << endl << i + 1 << "-ю строку делим на " << U[i][i] << endl;

		// Диагональный элемент делаем 1
		U[i] = U[i] / U[i][i];

		fout << U;

		// Занулляем всё, что под диагональю 
		for (int j = i + 1; j < N; j++)
		{
			fout << "\nИз " << j + 1 << "-й строки вычитаем " << i + 1 << "-ю, умноженную на " << U[j][i] << endl;

			U[j] = U[j] - U[i] * U[j][i];
			
			fout << U;
		}
	}
	return U;
}

// Создаём матрицу L
myMatrix CreateMatrixL(myMatrix A, myMatrix U, ofstream &fout)
{
	auto L = EmptyMatrix();
	fout << "\nВычисляем матрицу L: \n\n";
	for (int i = 0; i < N; i++)
	{
		fout << "ШАГ " << i + 1 << endl << endl;
		for (int j = 0; j <= i; j++)
		{
			fout << "Ставим на позицию (" << i + 1 << "," << j + 1 << ") матрицы L, элемент " << A[i][j] << ", \nкоторый берем с той же позиции матрицы А \n";
			L[i][j] = A[i][j]; 
			fout << L << endl;
			for (int k = 0; k < j; k++)
			{
				fout << "Вычитаем из L("<< i + 1 << "," << j + 1 << ") = "<< L[i][j] << " произведение L(" << i + 1 << "," << k + 1 << ")*U(" << k + 1 << "," << j + 1 << ") = (" << L[i][k] <<")*("<< U[k][j] << ")\n";
				L[i][j] -= L[i][k] * U[k][j];
				fout << L << endl;
			}
		}
	}
	return L;
}

// Вычисление определителя 
double Determinant(myMatrix L)
{
	double det = 1;
	for (int i = 0; i < N; i++)
		det *= L[i][i];
	return det;
}

// Создаем матрицу перестановок P
myMatrix CreateMatrixP()
{
	auto P = EmptyMatrix();
	for (int i = 0; i < N; i++)
		P[i][i] = 1;
	return P;
}

// Решаем систему, находим вектор X
vector<double> SolveX(vector<double> &b, myMatrix P, myMatrix L, myMatrix U)
{
	// Делаем перестановку элеметов вектора b в соответствии с матрицей перестановок
	auto b1 = P * b;
	
	// Решение уравнение Ly=b (Находим вектор y)
	vector<double> y; 
 	for (int i = 0; i < N; i++)
	{
		// Вычитаем из правой части все уже вычесленные значения y, домноженнные на коэффициенты
		for (int j = 0; j < i; j++)
			b1[i] -= L[i][j] * y[j];

		// Делим полученную правую часть на коэффициент при даигонали и ккладуём в вектор y
		y.push_back(b1[i] / L[i][i]);
	}

	// Решение уравнение Ux=y (находим вектор x)
	// Так как на диагоналях стоят 1, то делить на коэффициент не нужно, поэтому присвоим 
	auto x = y;				
	for (int i = N - 1; i >= 0; i--)
		for (int j = N - 1; j > i; j--)
			x[i] -= U[i][j] * x[j];
	return x;
}

// Находим вектор правых частей 
vector<double> SolveB(vector<double> x, myMatrix A)
{
	vector<double> b(N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			b[i] += A[i][j] * x[j];
	return b;
}

// Транспонирование матрицы
myMatrix TransposeMatrix(myMatrix M)
{
	myMatrix res = M;
	for (int i = 0; i < N; i++)
		for (int j = i + 1; j < N; j++)
		{
			double tmp = res[i][j];
			res[i][j] = res[j][i];
			res[j][i] = tmp;
		}
	return res;
}

// Создание обратной матрицы
myMatrix ReverseMatrix(myMatrix P, myMatrix L, myMatrix U, ofstream &fout) 
{
	fout << "Вычисление обратной матрицы: \n";
	myMatrix A_ = EmptyMatrix();
	// Единичнвя матрица
	myMatrix E = CreateMatrixP();
	
	fout << "Единичная матрица: \n" << E << endl; 
	
	for (int i = 0; i < N; i++)
	{
		fout << "ШАГ " << i + 1 << "\n\nРешим исходную систему, где вектор правых частей равен " << i+1 << "-й строке из матрицы E\nЗапишем полученный вектор в строку результирующей матрицы\n"; 
		// Решение исходной ситемы уравнений, где вектор b = строке из Е 
		A_[i] = SolveX(E[i], P, L, U);

		fout << A_ << endl;
	}

	fout << "ШАГ " << N + 1 << ": транспонируем полученную матрицу\n";
	return TransposeMatrix(A_);
}

// Максимальная сумма столбцов/строк в матрице
double MaxSumRowOrColumn(myMatrix A, int numNorm)
{
double max(0), sum(0);
for (int i = 0; i < A.size(); i++)
{
for (int j = 0; j < A.size(); j++)
	if (numNorm == 1)
		sum += fabs(A[i][j]);
	else 
		sum += fabs(A[j][i]);
if (sum > max)
	max = sum;
sum = 0;
}
return max;
}

// Вычисляем норму I
double NormI(myMatrix A)
{
return MaxSumRowOrColumn(A, 1);
}

// Вычисляем норму II
double NormII(myMatrix A)
{
return MaxSumRowOrColumn(A, 2);
}

// Вычисляем норму III
double NormIII(myMatrix A)
{
	// Эрмитово сопряженная матрица А*. Для действ. чисел просто транспонированная 
	auto At = TransposeMatrix(A);

	// А* x А
	auto AtA = At * A;

	real_2d_array a, vl, vr;
	a.setlength(N, N);

	real_1d_array  wl, wi;
	wl.setlength(N);
	wi.setlength(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = AtA[j][i];
	rmatrixevd(a , N, 0, wl, wi, vl, vr);

	double maxVal(DBL_MIN); 

	for (int i = 0; i < N; i++)
		if (fabs(wl[i]) > maxVal)
			maxVal = fabs(wl[i]);
	
	return maxVal;
}

// Основная ф-я 
int main()
{
	setlocale(LC_ALL, "Russian");

	// Файл для записи
	ofstream fout;
	
	// Файл для считывания 
	ifstream fin;

	// A - исх. матрица, AnotP - исходная матрица, строки которой не переставляем, LU = L*U, PA = P*A, _A - Обратная матрица
	myMatrix A, AnotP, L, U, LU, PA, PA_LU, A_, AA_, AA__E, P(CreateMatrixP()), E(CreateMatrixP());

	// Исходный вектор b и вектор решений х  
	vector<double> b, x, Ax, Ax_b;

	// Ранг матрицы A (предпологаем, что он равен 4-м)
	int rang = N; 

	// Определитель матрицы A
	double k = 1;

	// Открываемм файл для считывания данных 
	fin.open("var4.txt");

	// Считываем матрицу целиком
	fin >> A;
	AnotP = A; 

	// Считываем вектор целиком
	fin >> x;

	// Закрываем файл для считывания
	fin.close();

	fout.open("output.txt");

	// Выведем матрицу и вектор
	fout << "Вектор x: \n" << x << endl;
	fout << "Исходная матрица А: \n" << A << endl;
	
	b = SolveB(x, A);

	fout << "Вектор b: \n" << b << endl;

	// Создадим матрицу перестановок
	P = CreateMatrixP();

	// Создадим матрицу U
	U = CreateMatrixU(A, P, fout, rang, k);

	if (rang != N)
	{
		fout << "Ранг матрицы А = " << rang << "\n\nПродолжать решение невозможно :(";
		return 0;
	}

	// Создадим матрицу L
	L = CreateMatrixL(A, U, fout);

	fout << "\nРанг исходной матрицы А: " << rang;

	fout << "\nОпределитель матрицы A: " << Determinant(L) * k;

	// Рекшим СЛАУ 
	x = SolveX(b, P, L, U);

	fout << "\n\nВектор x: \n" << x << endl;

	A_ = ReverseMatrix(P, L, U, fout);
	fout << endl << A_ << endl;

	fout << "Итак, имеем: \n=============== \nL:\n" << L << "\nU:\n" << U << "\nP:\n" << P << "\nA:\n" << A << "\nA^(-1):\n" << A_ << "\nb:\n" << b << "\nx:\n" << x << "===============\n\n";

	PA = P * AnotP;
	fout << "P*A:" << endl << PA << endl;

	LU = L * U;
	fout << "L*U:" << endl << LU << endl;

	PA_LU = PA - LU;
	fout << "PA - LU:" << endl << PA_LU << "\nНорма несвязки: ||PA - LU|| = " << NormII(PA_LU) << endl << endl;

	Ax = AnotP * x;
	Ax_b = Ax - b;
	fout << "Ax: " << Ax << "\nAx-b: " << Ax_b;

	// Найдём норму вектора:
	double sum(0);
	for (int i = 0; i < N; i++)
		sum += fabs(Ax_b[i]);
	fout << "\nНорма несвязки: ||Ax - b|| = " << sum;

	AA_ = AnotP * A_;
	fout << "\n\nA*A^(-1):\n" << AA_ << endl;

	AA__E = AA_ - E;
	fout << "A*A^(-1) - E:\n" << AA__E << "\nНорма несвязки: ||A*A^(-1) - E|| = " << NormII(AA__E) << "\n\n";

	fout << "Число обусловленности исходной матррицы \n";
	fout << setw(25) << setiosflags(ios::left) << "В кубической норме: " << setprecision(6)  << NormI(A) * NormI(A_) << endl;
	fout << setw(25) << "В октаэдрической норме: "  << NormII(A) * NormII(A_) << endl;
	fout << setw(25) << "В евклидовой норме: " << sqrt(NormIII(A)) * sqrt(NormIII(A_)) << endl;

	fout.close();

	system("pause");
	return 0;
}