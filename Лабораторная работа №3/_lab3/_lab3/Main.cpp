#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <functional>
#include "linalg.h"

using namespace std;

int const N = 4;
double const EPS = 0.0001;
typedef vector<vector<double>> myMatrix;

#pragma region OldProgram

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
	//out << endl;
	return out;
}

// Умножение вектора на вектор (переопределение оператора *)
double operator * (vector<double> v1, vector<double> v2)
{
	double sum(0);
	for (int i = 0; i < N; i++)
		sum += v1[i] * v2[i];
	return sum;
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

// Сложение 2-х векторов поэлементно (переопределение оператора +)
vector<double> operator	+ (vector<double> v1, vector<double> v2)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v1[i] + v2[i]);
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
		out << fixed << setprecision(7) << setw(13) << v[i] << endl;
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
	vector<double> res(N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			res[i] += M[i][j]*v[j];
	return res;
}

// Перегрузка оператора - для вычитания двух матриц
myMatrix operator - (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			res[i][j] += M1[i][j] - M2[i][j
			];
	return res;
}

// Перегрузка оператора + для сложения двух матриц
myMatrix operator + (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			res[i][j] += M1[i][j] + M2[i][j];
	return res;
}

// Перегрузка оператора * для умножения матрицы на число
myMatrix operator * (myMatrix &M1, double &d)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			res[i][j] += M1[i][j] * d;
			
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
myMatrix CreateMatrixU(myMatrix &A, myMatrix &P, int &rang)
{
	auto U = A;
	int ind;

	// Пробегаем по всем строкам
	for (int i = 0; i < N; i++)
	{
		// Ищем индекс максимального элемента в столбце
		ind = IndexMaxNumInColumn(U, i);

		// Значит весь столбец (начиная со строки i и до N) состоит из 0, а значит ранг равен пройденному числу строк
		if (ind == -1)
		{
			rang = i;
			return U;    // Что-то вернуто-то надо
		}

		// Если найденный индекс не явл. индексом тек. строки, то меняем строки
		if (ind != i)
		{
			Swap(U, i, ind);
			Swap(A, i, ind);
			Swap(P, i, ind);
		}

		// Диагональный элемент делаем 1
		U[i] = U[i] / U[i][i];

		// Занулляем всё, что под диагональю 
		for (int j = i + 1; j < N; j++)
			U[j] = U[j] - U[i] * U[j][i];
	}
	return U;
}

// Создаём матрицу L
myMatrix CreateMatrixL(myMatrix A, myMatrix U)
{
	auto L = EmptyMatrix();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			L[i][j] = A[i][j];
			for (int k = 0; k < j; k++)
				L[i][j] -= L[i][k] * U[k][j];
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
myMatrix ReverseMatrix(myMatrix P, myMatrix L, myMatrix U)
{
	myMatrix A_ = EmptyMatrix();
	// Единичнвя матрица
	myMatrix E = CreateMatrixP();

	for (int i = 0; i < N; i++)
		A_[i] = SolveX(E[i], P, L, U);

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

// норма матрицы I
double NormI(myMatrix A)
{
	return MaxSumRowOrColumn(A, 1);
}

// норма матрицы II
double NormII(myMatrix A)
{
	return MaxSumRowOrColumn(A, 2);
}

// норма матрицы III
double NormIII(myMatrix A)
{
	auto At = TransposeMatrix(A);

	auto AtA = At * A;

	alglib::real_2d_array a, vl, vr;
	a.setlength(N, N);

	alglib::real_1d_array  wl, wi;
	wl.setlength(N);
	wi.setlength(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = AtA[j][i];
	rmatrixevd(a, N, 0, wl, wi, vl, vr);

	double maxVal(DBL_MIN);

	for (int i = 0; i < N; i++)
		if (fabs(wl[i]) > maxVal)
			maxVal = fabs(wl[i]);

	return maxVal;
}

// 1-я норма вектора
double NormVecI(vector <double> v)
{
	double m = 0;
	for (int i = 0; i < N; i++) {
		m = m < abs(v[i]) ? abs(v[i]) : m;
	}
	return m;
}

// 2-я норма вектора
double NormVecII(vector <double> v)
{
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += fabs(v[i]);
	}
	return sum;
}

// 3-я норма вектора
double NormVecIII(vector <double> v)
{
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += v[i] * v[i];
	}
	return sqrt(sum);
}

// Минимальная из норм вектора
double MinNormVec(vector<double> v)
{
	double N1 = NormVecI(v);
	double N2 = NormVecII(v);
	double N3 = NormVecIII(v);
	if (N1 < N2)
	{

		if (N1 < N3)
			return N1;
		else
			return N3;
	}
	else
		return N2;
}

#pragma endregion


// Метод простой итерации
vector<double> SimpleIterationMethod(myMatrix A, vector<double> b, int &factCountIter, ofstream &fout, unsigned int &time)
{
	unsigned int start_time = clock();

	fout << "МЕТОД ПРОСТОЙ ИТЕРАЦИИ\n\n";

	// xPrev соответствует Х(к-1), xNext - Х(к-1), xCur - Х(к)
	vector <double> xPrev(N), xCur(N), xNext(N), ax;

	double q(0), t(1.8 / NormII(A)), normDiscr;		// discrepancy - невязка
	// Счётчик итераций
	int countIter(0);
	fout << "Значения итерационного параметра: " << setprecision(7) << t << endl;
	fout << "Начальное приближение: " << b << endl << endl;
	fout << setw(44) << "Норма  " << setw(13) << "Оценка\n";
	fout << setw(5)  << "№" << setw(13) << "tau   " << setw(13) <<  "q   " << setw(13)
		 <<"невязки"<< setw(13) << "погр-ти" << setw(13) << setw(13) << "X[1]  " << setw(13) << "X[2]  " << setw(13) << "X[3]  "
		 << setw(13) << "X[4]  " << "\n\n";

	while (true)
	{
		ax = A * xCur;
		countIter++;

		// Вычисляем х(к+1)
		xNext = xCur + (b - ax ) * t;

		// Или каждый элемент 
		// for (int i = 0; i < N; i++) xNext[i] = xCur[i] + t * (b[i] - ax[i]); 

		if (countIter != 1)
			q = MinNormVec(xNext - xCur) / MinNormVec(xCur - xPrev);
		else q = MinNormVec(xNext);

		normDiscr = MinNormVec(ax - b) / MinNormVec(xNext);

		fout << setw(5) << countIter << setw(13) << t << setw(13) << q << setw(13) << normDiscr 
			 << setw(13) << MinNormVec(xNext-xCur) << xNext << "\n";
																															 
		// Если норма невязки меньше 10^(-4) то прекращаем работу цикла
		if (normDiscr < EPS)
			break;

		// Меняем предыдущее на текущее и текущее на следующее
		xPrev = xCur;
		xCur = xNext;
	}

	fout << "\nПолучили: \nЧисло итераций: " << countIter 
		 << "\nНорма невязки: " << normDiscr << "\n" 
		 << "Оценка нормы матрицы перехода q: " << q 
		 << "\nПослднее приближение: " << xNext << endl;

	factCountIter = countIter;

	unsigned int end_time = clock();
	time = end_time - start_time;

	return xNext;
}


// Градиентный метод наискорейшего спуска
vector<double> GradientDescentMethod(myMatrix A, vector<double> b, int &factCountIter, ofstream &fout, unsigned int &time)
{
	unsigned int start_time = clock();

	fout << "\n\nГРАДИЕНТНЫЙ МЕТОД НАИСКОРЕЙШЕГО СПУСКА\n\n";
	vector <double> xPrev(N), xCur(N), xNext(N), ax, r, ar;
	double q(0), t, normDiscr;		
	int countIter(0);

	fout << "Начальное приближение: " << b << endl << endl;
	fout << setw(44) << "Норма  " << setw(13) << "Оценка\n";
	fout << setw(5) << "№" << setw(13) << "tau   " << setw(13) << "q   " << setw(13)
		<< "невязки" << setw(13) << "погр-ти" << setw(13) << setw(13) << "X[1]  " << setw(13) << "X[2]  " << setw(13) << "X[3]  "
		<< setw(13) << "X[4]  " << "\n\n";

	while (true)
	{
		countIter++;
		ax = A * xCur;
		r = ax - b;
		ar = A * r;
		t = (r * r) / (ar * r);
		xNext = xCur - r * t;

		if (countIter != 1)
			q = MinNormVec(xNext - xCur) / MinNormVec(xCur - xPrev);
		else q = MinNormVec(xNext);

		normDiscr = MinNormVec(ax - b) / MinNormVec(xNext);


		fout << setw(5) << countIter << setw(13) << t << setw(13) << q << setw(13) << normDiscr
			 << setw(13) << MinNormVec(xNext - xCur) << xNext << "\n";

		if (normDiscr < EPS)
			break;

		xPrev = xCur;
		xCur = xNext;
	}

	fout << "Получили: \nЧисло итераций: " << countIter << "\nНорма невязки: " << normDiscr << "\n" << "Оценка нормы матрицы перехода q: " << q << "\nПослднее приближение: " << xNext << endl;

	factCountIter = countIter;

	unsigned int end_time = clock();
	time = end_time - start_time;

	return xNext;
}

// Вычсление следующей итерации для метода ПВР
vector <double> CalcNextX(myMatrix A, vector<double> b, vector<double> xCur, double w)
{
	vector<double> xNext(N);
	double sum;
	for (int i = 0; i < N; i++)
	{
		sum = 0;
		for (int j = 0; j < i; j++)
			sum += A[i][j] * xNext[j];

		for (int j = i + 1; j < N; j++)
			sum += A[i][j] * xCur[j];

		xNext[i] = xCur[i] + w*((b[i] - sum) / A[i][i] - xCur[i]);
	}
	return xNext;
}

// Поиск оптимального значения w
double wSearch(myMatrix A, vector<double> b, ofstream &fout)
{
	double const eps = 0.01;
	double w(0.1), wOptimal(0.1), s;
	int countIter, minCountIter(INT_MAX);
	for (w; w < 2; w += 0.1)
	{
		countIter = 0;
		vector<double> xCur(N), xNext(N);
		do
		{
			xCur = xNext;
			countIter++;
			xNext = CalcNextX(A, b, xCur, w);

		} while (MinNormVec(A*xCur - b) / MinNormVec(xNext)>= eps);

		if (countIter < minCountIter)
		{
			wOptimal = w;
			minCountIter = countIter;
		} 
		fout << "w = " << setprecision(1) << w << " - число итераций: " << countIter << endl;
	}
	fout << "\nМинимальное число итераций достигается при w = " << wOptimal << endl << endl;
	return wOptimal;
}

// Метол ПВР
vector<double> MethodSOR(myMatrix A, vector<double> b, int &factCountIter, ofstream &fout, unsigned int &time)
{
	unsigned int start_time = clock();

	fout << "\n\nМЕТОД ПВР(SOR)\n\n";
	vector <double> xPrev(N), xCur(N), xNext(N), ax;
	double q(0), s, normDiscr, w = wSearch(A, b, fout);		// discrepancy - невязка
	int countIter(0);

	fout << "Начальное приближение: " << b << endl << endl;
	
	fout << setw(44) << "Норма  " << setw(13) << "Оценка\n";
	fout << setw(5) << "№" << setw(13) << "w   " << setw(13) << "q   " << setw(13)
		<< "невязки" << setw(13) << "погр-ти" << setw(13) << setw(13) << "X[1]  " << setw(13) << "X[2]  " << setw(13) << "X[3]  "
		<< setw(13) << "X[4]  " << "\n\n";

	while (true)
	{
		countIter++;
		
		xNext = CalcNextX(A, b, xCur, w);

		ax = A * xCur;

		if (countIter > 1)
			q = MinNormVec(xNext - xCur) / MinNormVec(xCur - xPrev);
		else q = MinNormVec(xNext - xCur);

		normDiscr = MinNormVec(ax - b)/MinNormVec(xNext);


		fout << setw(5) << countIter << setw(13) << w << setw(13) << q << setw(13) << normDiscr
			 << setw(13) << MinNormVec(xNext - xCur) << xNext << "\n";

		if (normDiscr < EPS)
			break;

		xPrev = xCur;
		xCur = xNext;
	}

	fout << "Получили: \nЧисло итераций: " << countIter << "\nНорма невязки: " << normDiscr << "\n" << "Оценка нормы матрицы перехода q: " << q << "\nПослднее приближение: " << xNext << endl;

	factCountIter = countIter;

	unsigned int end_time = clock();
	time = end_time - start_time;

	return xNext;
}

// Метод сопряженных градиентов
vector<double> ConjugateGradientMethod(myMatrix A, vector<double> b,  int &factCountIter, ofstream &fout, unsigned int &time)
{
	unsigned int start_time = clock();

	fout << "\n\nМЕТОД СОПРЯЖЕННЫХ ГРАДИЕНТОВ\n\n";
	vector <double> xPrev(N), xCur(N), xNext(N), ax, r, rPrev, ar;
	double q(0), tCur, tNext(0), normDiscr, a(1);		// discrepancy - невязка
	int countIter(0);

	fout << "Начальное приближение: " << b << endl << endl;
	fout << setw(44) << "Норма  " << setw(13) << "Оценка\n";
	fout << setw(5) << "№" << setw(13) << "tau   " << setw(13) << "q   " << setw(13)
		<< "невязки" << setw(13) << "погр-ти" << setw(13) << setw(13) << "X[1]  " << setw(13) << "X[2]  " << setw(13) << "X[3]  "
		<< setw(13) << "X[4]  " << "\n\n";

	while (true)
	{
		countIter++;
		ax = A * xCur;
		r = ax - b; //вычисляем вектор невязки
		tCur = tNext; //переопред t для вычисления альфа
		ar = A * r;
		tNext = (r * r) / (ar * r);
		if (countIter > 1)
		{
			a = 1 / (1 - (tNext / (tCur*a))*((r * r) / (rPrev * rPrev)));
		}
		xNext = xCur * a + xPrev * (1 - a) - r * tNext*a;

		if (countIter != 1)
			q = MinNormVec(xNext - xCur) / MinNormVec(xCur - xPrev);
		else q = MinNormVec(xNext);

		normDiscr = MinNormVec(A * xCur - b) / MinNormVec(xNext);

		fout << setw(5) << countIter << setw(13) << tNext << setw(13) << q << setw(13) << normDiscr
		     << setw(13) << MinNormVec(xNext - xCur) << xNext <<  "\n";

		rPrev = r;
		xPrev = xCur;
		xCur = xNext;

		if (normDiscr < EPS)
			break;
	}

	fout << "Получили: \nЧисло итераций: " << countIter << "\nНорма невязки: " << normDiscr << "\n" << "Оценка нормы матрицы перехода q: " << q << "\nПослднее приближение: " << xNext << endl;

	factCountIter = countIter;

	unsigned int end_time = clock();
	time = end_time - start_time;

	return xNext;
}

// Основная ф-я 
int main()
{
	setlocale(LC_ALL, "Russian");

	// Номер набора
	int num = 14;

	// Файл для записи
	ofstream fout;

	// Файл для считывания 
	ifstream fin;

	// A - исх. матрица
	myMatrix A, L, U, P(CreateMatrixP());

	// Исходный вектор b и вектор решений х  
	vector<double> b, x;

	// Открываемм файл для считывания данных 
	fin.open("input" + to_string(num) + ".txt");

	// Считываем матрицу целиком
	fin >> A;

	// Считываем вектор целиком
	fin >> b;

	// Ранг матрицы A (предпологаем, что он равен 4-м)
	int rang = N;

	// Закрываем файл для считывания
	fin.close();

	//fout.open("output.txt");
	fout.open("output" + to_string(num) + ".txt");

	// Выведем матрицу и вектор
	fout << "Вектор b: \n" << b << endl << endl;
	fout << "Исходная матрица А: \n" << A << endl;

	// Норма матрицы
	fout << "Норма матрицы А: " << NormII(A) << endl << endl;

	// Фактическое число итераций в каждом из методов
	int mFactSIM, mFactGDM, mFactSOR, mFactCGM;

	// Время работы каждого из методов
	unsigned int timeSIM, timeGDM, timeSOR, timeCGM;

	// Метод простой итерации
	auto xSIM = SimpleIterationMethod(A, b, mFactSIM, fout, timeSIM);

	// Градиентный метод наискорейшего спуска
	auto xGDM = GradientDescentMethod(A, b, mFactGDM, fout, timeGDM);

	//Метод ПВР
	auto xSOR = MethodSOR(A, b, mFactSOR, fout, timeSOR);

	// Метод сопряженных градиентов
	auto xCGM = ConjugateGradientMethod(A, b, mFactCGM, fout, timeCGM);

	// LU разложение
	U = CreateMatrixU(A, P, rang);
	L = CreateMatrixL(A, U);

	// Метод LU разложения 
	auto xLU = SolveX(b, P, L, U);

	// Число обусловленности (для сравнения числа итераций)
	double cond(NormII(A)*NormII(ReverseMatrix(P, L, U)));


	//Теоретическое исло итераций по каждому из методов
	int mTeorSIM(round(0.5 * cond * log(1 / EPS))),
		mTeorGDM(round(0.5 * cond * log(1 / EPS))),
		mTeorSOR(round(0.25 * sqrt(cond) * log(1 / EPS))),
		mTeorCGM(round(sqrt(cond) * log(2 / EPS) * 0.5));

	fout << "\n\nПОДВЕДЕМ ИТОГИ: \n\n";
	fout << "Решение методом LU разложения: \n" << xLU << endl << endl;

	fout << "\nЧисло обусловленности: " << cond << endl;

	fout << "Решение методом простой итерации: \n" << xSIM << endl << 
			"Фактическое число итераций:   " << mFactSIM << 
			"\nТеоретическое число итераций: " << mTeorSIM << 
		    "\nСкорость работы метода: " << timeSIM / 1000.0 << " сек.\n\n";

	fout << "Решение градиентным методом наискорейшего спуска: \n" << xGDM << endl <<
			"Фактическое число итераций:   " << mFactGDM <<
			"\nТеоретическое число итераций: " << mTeorGDM <<
			"\nСкорость работы метода: " << timeGDM / 1000.0 << " сек.\n\n";

	fout << "Решение методом ПВР: \n" << xSOR << endl <<
			"Фактическое число итераций:   " << mFactSOR <<
			"\nТеоретическое число итераций: " << mTeorSOR <<
			"\nСкорость работы метода: " << timeSOR / 1000.0 << " сек.\n\n";

	fout << "Решение методом сопряженных градиентов: \n" << xCGM << endl <<
			"Фактическое число итераций:   " << mFactCGM <<
			"\nТеоретическое число итераций: " << mTeorCGM <<
			"\nСкорость работы метода: " << timeCGM / 1000.0 << " сек.\n\n";


	fout.close();

	//system("pause");
	return 0;
}