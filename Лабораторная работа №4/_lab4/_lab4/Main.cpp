#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <functional>

using namespace std;

double const EPS = 0.0001;

struct Equation
{
	double multiplier;
	double multiplier_1;   // Коэффициент при выражаемой переменной (х или у)
	char variable;		 // Выражаемая переменная
	string trigFunc;	 // Тригономметрическая функция 
	char variableInTrigFunk;
	double addInTrigFunk;
	char signForTrig = '+';
	double freeMember;   // Свободный член (в правой части)
	string f;			 // Вид функции f
} eq1, eq2;

void Parser(string s, Equation &eq)
{
	int i = 0;
	string f = "";
	while (s[i] != '\0')
	{
		if (s[i] == '=')
		{
			i++;
			if (s[i] == '-')
			{
				f += '+';
				i++;
			}
			else f += '-';
		}
		f += s[i];
		i++;
	}
	f += "=0";
	eq.f = f;
	i = 0;
	while (s[i] != '\0')
	{
		if (s[i] == 's' || s[i] == 'c')
		{
			if (i == 0)
			{
				eq.signForTrig = '-';
			}
			else if (s[i - 1] == '-')
			{
				eq.signForTrig = '+';
			}
			else eq.signForTrig = '-';

			eq.trigFunc = s.substr(i, 3);
			i+=4;
			eq.variableInTrigFunk = s[i];
			i++;
			if (s[i] == ')')
			{
				eq.addInTrigFunk = 0;
				i++;
			}
			else if (s[i+2] == ')')
			{
				eq.addInTrigFunk = stof(s.substr(i, 2));
				i += 3;
			}
			else
			{
				cout <<s.substr(i, 4);
				eq.addInTrigFunk = stof(s.substr(i, 4));
				i += 5;
			}
		}
		else if (s[i] <= '9' && s[i] >= '0')
		{
			if (s[i + 1] == 'x' || s[i + 1] == 'y')
			{
				if (i != 0 && s[i - 1] == '-')
					eq.multiplier = stof(s.substr(i - 1, 2));
				else eq.multiplier = stof(s.substr(i, 1));
				i ++;
			}
			else
			{
				string freeMember = "";
				while (s[i] != '\0')
				{
					freeMember += s[i];
					i++;
				}
				eq.freeMember = stof(freeMember);
			}
		}
		else if (s[i] == 'x' || s[i] == 'y')
		{
			eq.variable = s[i];
			if (i != 0)
			{
				if (s[i - 1] == '-')
					eq.multiplier = -1;
				else if (s[i - 1] == '+')
					eq.multiplier = 1;
			}
			else eq.multiplier = 1;
			i++;
		}
		else i++;
	}
	eq.multiplier_1 = 1 / eq.multiplier;
}
 
double NormVecI(double x, double y)
{
	return max(fabs(x), fabs(y));
}
// 2-я норма вектора
double NormVecII(double x, double y)
{
	return fabs(x)+fabs(y);
}

// 3-я норма вектора
double NormVecIII(double x, double y)
{
	return sqrt(x*x+y*y);
}

double SolveXY(Equation eq, double var)
{
	double res = eq.freeMember;
	double trig = eq.addInTrigFunk + var;
	if (eq.trigFunc[0] == 's')
		trig = sin(trig);
	else trig = cos(trig);
	if (eq.signForTrig == '-')
		trig = 0 - trig;
	res += trig;
	res *= eq.multiplier_1;

	return res;
}

double fSolve(Equation eq, double var1, double var2)
{
	double res = eq.multiplier*var2;
	double trig = eq.addInTrigFunk + var1;
	if (eq.trigFunc[0] == 's')
		trig = sin(trig);
	else trig = cos(trig);
	if (eq.signForTrig == '+')
		trig = 0 - trig;
	res += trig;
	res -= eq.freeMember;

	return res;
}

double SolveJacobian(Equation eq, double var)
{
	double trig = eq.addInTrigFunk + var;
	if (eq.trigFunc[0] == 's')
		trig = cos(trig);
	else trig = 0-sin(trig);
	if (eq.signForTrig == '-')
		trig = 0 - trig;
	return eq.multiplier_1*trig;
}

void PrintJacobian(Equation eq, ofstream & fout)
{
	if (abs(eq.multiplier_1) != 1) fout << eq.multiplier_1 << "*";
	if (eq.trigFunc == "cos")
	{
		if (eq.signForTrig == '-' && eq.multiplier_1 < 0 || eq.signForTrig == '+' && eq.multiplier_1 > 0)
			fout << "-sin(" << eq.variableInTrigFunk;
		else fout << "sin(" << eq.variableInTrigFunk;
	}
	else
	{
		if (eq.signForTrig == '-' && eq.multiplier_1 < 0 || eq.signForTrig == '+' && eq.multiplier_1 > 0)
			fout << "cos(" << eq.variableInTrigFunk;
		else fout << "-cos(" << eq.variableInTrigFunk;
	}
	if (eq.addInTrigFunk != 0)
		if (eq.addInTrigFunk < 0)
			fout << eq.addInTrigFunk;
		else fout << "+" << eq.addInTrigFunk;
	fout << ')';
}

void Jacobian(ofstream & fout)
{
	fout << "Якобиан: \n\n";
	PrintJacobian(eq1, fout);
	fout << setw(10) << 0 << endl << 0 << setw(15);
	PrintJacobian(eq2, fout);
}

void PrintPhi(Equation eq, ofstream &fout)
{
	fout << eq.variable << "=";
	if (eq.multiplier_1 != 1)
		fout << eq.multiplier_1 << "*" << "(";
	fout << eq.freeMember << eq.signForTrig << eq.trigFunc << "(" << eq.variableInTrigFunk;
	if (eq.addInTrigFunk != 0)
		if (eq.addInTrigFunk < 0)
			fout << eq.addInTrigFunk;
		else fout << "+" << eq.addInTrigFunk;
	fout << ")";
	if (eq.multiplier_1 != 1)
		fout << ")";
}

tuple<double, double> SimpleIteration(ofstream & fout, double x0, double y0)
{
	fout << "\n\nМЕТОД ПРОСТОЙ ИТЕРАЦИИ: \n";
	fout << "Представим уравнения системы в виде f(x,y) = 0: \n";
	fout << "\nУр. №1: " << eq1.f << "\nУр. №2: " << eq2.f << endl;
	fout << "\n\nПреобразуем систему к эквивалентной:";
	fout << "\nУр. №1: ";
	PrintPhi(eq1, fout);
	fout << "\nУр. №2: ";
	PrintPhi(eq2, fout);

	double x, y, xNext(x0), yNext(y0), f1, f2;
	int countItr = 0;
	fout << "\n\nНачальное приближение: x = " <<  x0  << ", y = " << y0 << "\n";
	Jacobian(fout);
	fout << "\nПодставим начальное приближение в якобиан: \n\n";
	double x_proizv = SolveJacobian(eq2, y0), y_proizv = SolveJacobian(eq1, x0);
	fout << setw(10) << x_proizv << setw(10) << 0 << endl << setw(10) << 0 << setw(10) << y_proizv;
	double Norm = max(abs(x_proizv), abs(y_proizv));
	fout << "\n\nНорма якобиана: " << Norm << endl << endl;
	do
	{
		countItr++;
		x = xNext;
		y = yNext;
		xNext = SolveXY(eq2, y);
		yNext = SolveXY(eq1, x);
		f1 = fSolve(eq1, xNext, yNext);
		f2 = fSolve(eq2, yNext, xNext);
		Norm = NormVecIII(x - xNext, y - yNext);
		fout << setw(3) << countItr << setw(13) << fixed << setprecision(7) << xNext << setw(13) << yNext << setw(13) 
			 << Norm << setw(20) <<  uppercase << scientific << f1 << setw(20) << f2;
		fout.copyfmt(ios(NULL));
		fout << setw(13) << NormVecI(SolveJacobian(eq2, yNext), SolveJacobian(eq1, xNext)) << "\n";
	} while (Norm >= EPS);

	fout << "Число итераций: " << countItr << endl;
	fout << "Последнее приближение: x = " << fixed << setprecision(7) << xNext << "  y = " << setw(13) << yNext << endl;
	return tuple<double, double>(xNext, yNext);
}

void df(Equation eq, ofstream & fout)
{
	if (eq.trigFunc[0] == 's')
	{
		if (eq.signForTrig == '+')
			fout << '-';
		fout << "cos(";
	}
	else
	{
		if (eq.signForTrig == '-')
			fout << '-';
		fout << "sin(";
	}
	fout << eq.variableInTrigFunk;
	if (eq.addInTrigFunk > 0)
		fout << '+';
	if (eq.addInTrigFunk != 0)
		fout << eq.addInTrigFunk;
	fout << ')';
}

double SolveDf(Equation eq, double var1)
{
	double d;
	if (eq.trigFunc[0] == 's') 
		d = cos(var1 + eq.addInTrigFunk);
	else
		d = 0 - sin(var1 + eq.addInTrigFunk);
	if (eq.signForTrig == '+')
		d = 0 - d;

	return d;
}

tuple<double, double> Newton(ofstream & fout, double x0, double y0)
{
	fout << "\n\nМЕТОД НЬЮТОНА: \n";
	fout << "Представим уравнения системы в виде f(x,y) = 0: \n";
	fout << "\nУр. №1: " << eq1.f << "\nУр. №2: " << eq2.f << endl;
	fout << "Матрица производных имеет вид: \n";
	fout.copyfmt(ios(NULL));
	df(eq1, fout);  fout << "\t" << eq1.multiplier << endl;
	fout << eq2.multiplier << "\t";
	df(eq2, fout); fout << endl;
	double x, y, xNext(x0), yNext(y0), f1, f2;
	double dxf1, dxf2, dyf1, dyf2, det, Norm;
	int countItr = 0;
	dxf2 = eq2.multiplier;
	dyf1 = eq1.multiplier;
	f1 = fSolve(eq1, x0, y0);
	f2 = fSolve(eq2, y0, x0);
	do 
	{
		countItr++;
		x = xNext;
		y = yNext;

		dxf1 = SolveDf(eq1, x);
		dyf2 = SolveDf(eq2, y);

		det = dxf1 * dyf2 - dxf2 * dyf1;
		xNext = x - (f1 * dyf2 - f2 * dyf1) / det;
		yNext = y - (f2 * dxf1 - f1 * dxf2) / det;
		
		f1 = fSolve(eq1, xNext, yNext);
		f2 = fSolve(eq2, yNext, xNext);

		Norm = NormVecIII(x - xNext, y - yNext);
		fout << setw(3) << countItr << setw(13) << fixed << setprecision(7) << xNext << setw(13) << yNext << setw(20)
			 << uppercase << scientific << Norm << setw(20)  << f1 << setw(20) << f2 << endl;
	} while (Norm >= EPS);

	fout << "Число итераций: " << countItr << endl;
	fout << "Последнее приближение: x = " << fixed << setprecision(7) << xNext << "  y = " << setw(13) << yNext << endl;
	return tuple<double, double>(xNext, yNext);
}


double F(double x, double y)
{
	double f1 = fSolve(eq1, x, y), f2 = fSolve(eq2, y, x);
	return f1 * f1 + f2 * f2;
}


double dFdx(double x, double y)
{
	return 2 * fSolve(eq1, x, y)*SolveDf(eq1, x) + 2 * fSolve(eq2, y, x)*eq2.multiplier;
}

double dFdy(double x, double y)
{
	return 2 * fSolve(eq1, x, y)*eq1.multiplier + 2 * fSolve(eq2, y, x)*SolveDf(eq2, y);
}


double GradIter(double x0, double y0)
{
	double xNext = x0, yNext = y0, a = 0, a0 = 1, l = 0.01, l1, l0 = 0.01, x = x0, y = y0;
	int iter = 0, min;
	min = INT_MAX;
	for (int i = 0; i < 80; i++)
	{
		xNext = x;
		yNext = y;
		iter = 0;
		do {
			iter++;
			x0 = xNext;
			y0 = yNext;
			a = a0;
			while (F(x0 - a * dFdx(x0, y0), y0 - a * dFdy(x0, y0)) >= F(x0, y0)) {
				a *= l;
			}
			xNext = x0 - a * dFdx(x0, y0);
			yNext = y0 - a * dFdy(x0, y0);

		} while (NormVecIII(x0 - xNext, y0 - yNext) >= EPS);
		if (iter < min)
		{
			l1 = l;
			min = iter;
		}
		l += l0;
	}
	return l1;
}

tuple<double, double> Gradient(ofstream &fout, double x0, double y0)//Градиентный спуск
{
	using namespace std::placeholders;
	double xNext(x0), yNext(y0), a, a0(1), l, k;
	l = GradIter(x0, y0);
	fout << "\nОптимальный параметр lambda = " << l << endl;
	int countItr = 0;
	do {
		countItr++;
		x0 = xNext;
		y0 = yNext;
		a = 0.5;
		k = 0;
		while (F(x0 - a * dFdx(x0, y0), y0 - a * dFdy(x0, y0)) >= F(x0, y0))
		{
			a *= l;
			k++;
		}
	
		double f1 = fSolve(eq1, x0, y0), f2 = fSolve(eq2, y0, x0);
		xNext = x0 - a * dFdx(x0, y0);
		yNext = y0 - a * dFdy(x0, y0);
		fout << setw(3) << countItr << setw(13) << fixed << setprecision(7) << xNext << setw(13) << yNext << setw(13) << a << setw(20)
			 << uppercase << scientific << NormVecIII(f1, f2) << setw(20) << f1 << setw(20) << f2 << setw(20) << NormVecIII(f1, f2)*NormVecIII(f1, f2) << setw(3) << k << endl;
	} while (NormVecIII(x0 - xNext, y0 - yNext) >= EPS);
	fout << "Число итераций: " << countItr << endl;
	fout << "Последнее приближение: x = " << fixed << setprecision(7) << xNext << "  y = " << setw(13) << yNext << endl;
	return tuple<double, double>(xNext, yNext);
}




int main()
{
	setlocale(LC_ALL, "Russian");

	// Номер набора
	int num = 9;

	// Начальное приближение
	double x0, y0;

	// Файл для считывания 
	ifstream fin;

	// Файл для записи
	ofstream fout;

	string strEq1, strEq2;

	// Открываемм файл для считывания данных 
	fin.open("input" + to_string(num) + ".txt");

	fin >> strEq1 >> strEq2;
	fin >> x0 >> y0;
	fout.open("output"+to_string(num) + ".txt");
	Parser(strEq1, eq1);
	Parser(strEq2, eq2);
	fout << "ИСХОДНАЯ СИСТЕМА УРАВНЕНИЙ: ";
	fout << "\nУр. №1: " << strEq1;
	fout << "\nУр. №2: " << strEq2;
	SimpleIteration(fout, x0, y0);
	Newton(fout, x0, y0);
	Gradient(fout, x0, y0);
	//system("pause");
	return 0;
}