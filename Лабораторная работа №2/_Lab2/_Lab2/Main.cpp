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

// ���������� �������� >> ��� �������
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

// ���������� �������� << ��� �������
ostream &operator << (ostream &out, vector<double> &v)
{
	for (int i = 0; i < N; i++)
		// 8 ������ ����� ������� � ������ ������� 15 ��������
		out << fixed << setprecision(7) << setw(13) << v[i];
	out << endl;
	return out;
}

// ��������� ������� �� ����� (��������������� ��������� *)
vector<double> operator * (vector<double> v, double num)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v[i] * num);
	return res;
}

// ������� ������� �� ����� (��������������� ��������� /)
vector<double> operator / (vector<double> v, double num)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v[i] / num);
	return res;
}

// �������� 2-� �������� ����������� (��������������� ��������� -)
vector<double> operator	- (vector<double> v1, vector<double> v2)
{
	vector<double> res;
	for (int i = 0; i < N; i++)
		res.push_back(v1[i] - v2[i]);
	return res;
}

// ������ ������ �������
myMatrix EmptyMatrix()
{
	vector<double> v(N);
	myMatrix emp;
	for (int i = 0; i < N; i++)
		emp.push_back(v);
	return emp;
}

// ���������� ��������� >> ��� �������
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

// ���������� ��������� << ��� �������
ostream &operator << (ostream &out, myMatrix &v)
{
	for (int i = 0; i < N; i++)
		// 8 ������ ����� ������� � ������ ������� 15 ��������
		out << fixed << setprecision(7) << setw(13) << v[i];
	return out;
}

// ���������� ��������� * ��� ������������ ���� ������
myMatrix operator * (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			for (int k = 0; k < M1[0].size(); k++)
				res[i][j] += M1[i][k] * M2[k][j];
	return res;
}

//  ���������� ��������� * ��� ��������� ������� �� ������ 
vector<double> operator * (myMatrix &M, vector<double> &v)
{
	vector<double> res(v.size());
	for (int i = 0; i < M.size(); i++)
		for (int j = 0; j < v.size(); j++)
			res[i] += v[j] * M[i][j];			
	return res;
}

// ���������� ��������� - ��� ������������ ���� ������
myMatrix operator - (myMatrix &M1, myMatrix &M2)
{
	myMatrix res = EmptyMatrix();
	for (int i = 0; i < M1.size(); i++)
		for (int j = 0; j < M1[0].size(); j++)
			res[i][j] += M1[i][j] - M2[i][j
];
	return res;
}

// ������� ������ ������� ������ � ������� ����� ������������ ������� ������� �������
int IndexMaxNumInColumn(myMatrix M, int columnIndex)
{
	int ind = columnIndex;
	for (int i = columnIndex; i < N; i++)
		if (fabs(M[i][columnIndex]) > fabs(M[ind][columnIndex]))
			ind = i;
	if (M[ind][columnIndex] == 0) // �������� 0 �� ���������
		return -1;
	return ind;
}

// �-�, ������� ������ ������� ������ 
void Swap(myMatrix &M, int index1, int index2)
{
	auto tmp = M[index1];
	M[index1] = M[index2];
	M[index2] = tmp;
}

// ������ ������� U
myMatrix CreateMatrixU(myMatrix &A, myMatrix &P, ofstream &fout, int &rang, double &k)
{
	auto U = A;
	fout << "���������� ������� U: " << endl;
	int ind;

	// ��������� �� ���� �������
	for (int i = 0; i < N; i++)
	{
		fout << "\n��� " << i + 1 << endl << endl;

		// ���� ������ ������������� �������� � �������
		ind = IndexMaxNumInColumn(U, i);

		fout << "������������ ������� � ������� " << i + 1 << ": " << fabs(U[ind][i]) << endl;

		// ������ ���� ������� (������� �� ������ i � �� N) ������� �� 0, � ������ ���� ����� ����������� ����� �����
		if (ind == -1)
		{
			rang = i;
			return U;    // ���-�� �������-�� ����
		}
		
		// ���� ��������� ������ �� ���. �������� ���. ������, �� ������ ������
		if (ind != i) 
		{
			k *= -1;
			Swap(U, i, ind);
			Swap(A, i, ind);
			Swap(P, i, ind);
			fout << "������ ������� " << i + 1 << "-� � " << ind + 1 << "-� ������ \n" << U;
		}
		else fout << "������������ ����� �� ��������� \n";

		fout << endl << i + 1 << "-� ������ ����� �� " << U[i][i] << endl;

		// ������������ ������� ������ 1
		U[i] = U[i] / U[i][i];

		fout << U;

		// ��������� ��, ��� ��� ���������� 
		for (int j = i + 1; j < N; j++)
		{
			fout << "\n�� " << j + 1 << "-� ������ �������� " << i + 1 << "-�, ���������� �� " << U[j][i] << endl;

			U[j] = U[j] - U[i] * U[j][i];
			
			fout << U;
		}
	}
	return U;
}

// ������ ������� L
myMatrix CreateMatrixL(myMatrix A, myMatrix U, ofstream &fout)
{
	auto L = EmptyMatrix();
	fout << "\n��������� ������� L: \n\n";
	for (int i = 0; i < N; i++)
	{
		fout << "��� " << i + 1 << endl << endl;
		for (int j = 0; j <= i; j++)
		{
			fout << "������ �� ������� (" << i + 1 << "," << j + 1 << ") ������� L, ������� " << A[i][j] << ", \n������� ����� � ��� �� ������� ������� � \n";
			L[i][j] = A[i][j]; 
			fout << L << endl;
			for (int k = 0; k < j; k++)
			{
				fout << "�������� �� L("<< i + 1 << "," << j + 1 << ") = "<< L[i][j] << " ������������ L(" << i + 1 << "," << k + 1 << ")*U(" << k + 1 << "," << j + 1 << ") = (" << L[i][k] <<")*("<< U[k][j] << ")\n";
				L[i][j] -= L[i][k] * U[k][j];
				fout << L << endl;
			}
		}
	}
	return L;
}

// ���������� ������������ 
double Determinant(myMatrix L)
{
	double det = 1;
	for (int i = 0; i < N; i++)
		det *= L[i][i];
	return det;
}

// ������� ������� ������������ P
myMatrix CreateMatrixP()
{
	auto P = EmptyMatrix();
	for (int i = 0; i < N; i++)
		P[i][i] = 1;
	return P;
}

// ������ �������, ������� ������ X
vector<double> SolveX(vector<double> &b, myMatrix P, myMatrix L, myMatrix U)
{
	// ������ ������������ �������� ������� b � ������������ � �������� ������������
	auto b1 = P * b;
	
	// ������� ��������� Ly=b (������� ������ y)
	vector<double> y; 
 	for (int i = 0; i < N; i++)
	{
		// �������� �� ������ ����� ��� ��� ����������� �������� y, ������������ �� ������������
		for (int j = 0; j < i; j++)
			b1[i] -= L[i][j] * y[j];

		// ����� ���������� ������ ����� �� ����������� ��� ��������� � ������� � ������ y
		y.push_back(b1[i] / L[i][i]);
	}

	// ������� ��������� Ux=y (������� ������ x)
	// ��� ��� �� ���������� ����� 1, �� ������ �� ����������� �� �����, ������� �������� 
	auto x = y;				
	for (int i = N - 1; i >= 0; i--)
		for (int j = N - 1; j > i; j--)
			x[i] -= U[i][j] * x[j];
	return x;
}

// ������� ������ ������ ������ 
vector<double> SolveB(vector<double> x, myMatrix A)
{
	vector<double> b(N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			b[i] += A[i][j] * x[j];
	return b;
}

// ���������������� �������
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

// �������� �������� �������
myMatrix ReverseMatrix(myMatrix P, myMatrix L, myMatrix U, ofstream &fout) 
{
	fout << "���������� �������� �������: \n";
	myMatrix A_ = EmptyMatrix();
	// ��������� �������
	myMatrix E = CreateMatrixP();
	
	fout << "��������� �������: \n" << E << endl; 
	
	for (int i = 0; i < N; i++)
	{
		fout << "��� " << i + 1 << "\n\n����� �������� �������, ��� ������ ������ ������ ����� " << i+1 << "-� ������ �� ������� E\n������� ���������� ������ � ������ �������������� �������\n"; 
		// ������� �������� ������ ���������, ��� ������ b = ������ �� � 
		A_[i] = SolveX(E[i], P, L, U);

		fout << A_ << endl;
	}

	fout << "��� " << N + 1 << ": ������������� ���������� �������\n";
	return TransposeMatrix(A_);
}

// ������������ ����� ��������/����� � �������
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

// ��������� ����� I
double NormI(myMatrix A)
{
return MaxSumRowOrColumn(A, 1);
}

// ��������� ����� II
double NormII(myMatrix A)
{
return MaxSumRowOrColumn(A, 2);
}

// ��������� ����� III
double NormIII(myMatrix A)
{
	// �������� ����������� ������� �*. ��� ������. ����� ������ ����������������� 
	auto At = TransposeMatrix(A);

	// �* x �
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

// �������� �-� 
int main()
{
	setlocale(LC_ALL, "Russian");

	// ���� ��� ������
	ofstream fout;
	
	// ���� ��� ���������� 
	ifstream fin;

	// A - ���. �������, AnotP - �������� �������, ������ ������� �� ������������, LU = L*U, PA = P*A, _A - �������� �������
	myMatrix A, AnotP, L, U, LU, PA, PA_LU, A_, AA_, AA__E, P(CreateMatrixP()), E(CreateMatrixP());

	// �������� ������ b � ������ ������� �  
	vector<double> b, x, Ax, Ax_b;

	// ���� ������� A (������������, ��� �� ����� 4-�)
	int rang = N; 

	// ������������ ������� A
	double k = 1;

	// ���������� ���� ��� ���������� ������ 
	fin.open("var4.txt");

	// ��������� ������� �������
	fin >> A;
	AnotP = A; 

	// ��������� ������ �������
	fin >> x;

	// ��������� ���� ��� ����������
	fin.close();

	fout.open("output.txt");

	// ������� ������� � ������
	fout << "������ x: \n" << x << endl;
	fout << "�������� ������� �: \n" << A << endl;
	
	b = SolveB(x, A);

	fout << "������ b: \n" << b << endl;

	// �������� ������� ������������
	P = CreateMatrixP();

	// �������� ������� U
	U = CreateMatrixU(A, P, fout, rang, k);

	if (rang != N)
	{
		fout << "���� ������� � = " << rang << "\n\n���������� ������� ���������� :(";
		return 0;
	}

	// �������� ������� L
	L = CreateMatrixL(A, U, fout);

	fout << "\n���� �������� ������� �: " << rang;

	fout << "\n������������ ������� A: " << Determinant(L) * k;

	// ������ ���� 
	x = SolveX(b, P, L, U);

	fout << "\n\n������ x: \n" << x << endl;

	A_ = ReverseMatrix(P, L, U, fout);
	fout << endl << A_ << endl;

	fout << "����, �����: \n=============== \nL:\n" << L << "\nU:\n" << U << "\nP:\n" << P << "\nA:\n" << A << "\nA^(-1):\n" << A_ << "\nb:\n" << b << "\nx:\n" << x << "===============\n\n";

	PA = P * AnotP;
	fout << "P*A:" << endl << PA << endl;

	LU = L * U;
	fout << "L*U:" << endl << LU << endl;

	PA_LU = PA - LU;
	fout << "PA - LU:" << endl << PA_LU << "\n����� ��������: ||PA - LU|| = " << NormII(PA_LU) << endl << endl;

	Ax = AnotP * x;
	Ax_b = Ax - b;
	fout << "Ax: " << Ax << "\nAx-b: " << Ax_b;

	// ����� ����� �������:
	double sum(0);
	for (int i = 0; i < N; i++)
		sum += fabs(Ax_b[i]);
	fout << "\n����� ��������: ||Ax - b|| = " << sum;

	AA_ = AnotP * A_;
	fout << "\n\nA*A^(-1):\n" << AA_ << endl;

	AA__E = AA_ - E;
	fout << "A*A^(-1) - E:\n" << AA__E << "\n����� ��������: ||A*A^(-1) - E|| = " << NormII(AA__E) << "\n\n";

	fout << "����� ��������������� �������� �������� \n";
	fout << setw(25) << setiosflags(ios::left) << "� ���������� �����: " << setprecision(6)  << NormI(A) * NormI(A_) << endl;
	fout << setw(25) << "� �������������� �����: "  << NormII(A) * NormII(A_) << endl;
	fout << setw(25) << "� ���������� �����: " << sqrt(NormIII(A)) * sqrt(NormIII(A_)) << endl;

	fout.close();

	system("pause");
	return 0;
}