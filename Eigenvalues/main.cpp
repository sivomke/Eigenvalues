#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

const int n = 10;
const double eps = 0.000000000000001;



double prod(double* x, double* y) {
	double prod = 0;
	for (int i = 0; i < n; ++i) {
		prod += x[i] * y[i];
	}
	return prod;
}

void matrix_max_ind(double(*A)[n], int* i, int* j) { //index of max off diagonal elem by abs value in matrix
	double max = abs(A[1][0]);
	for (int m = 0; m < n; ++m) {
		for (int k = 0; k < n; ++k) {
			if ((abs(A[m][k]) >= max) && (k != m)) {
				max = abs(A[m][k]);
				*i = m;
				*j = k;
			}
		}
	}

}

void print_matrix(double(*A)[n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << setw(15) << right << scientific << A[i][j];
		}
		std::cout << endl;
	}
	std::cout << endl;
}

double sgn(double x) {
	if (x > 0) return 1.;
	else return -1.;
}

double check_conv(double(*A)[n]) {
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j != i) {
				sum += A[i][j] * A[i][j];  //change to squares!!! sum += A[i][j] * A[i][j]
			}
		}
	}
	return sum;
}

double Jacoby_eigenv(double(*A)[n]) {
	double A_cur[n][n] = { 0 };
	double B[n][n] = { 0 };
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B[i][j] = A[i][j];
		}
	}
	double T, t, c, s, error;
	int* a = new int(0);
	int* b = new int(0);
	do {
		matrix_max_ind(B, a, b);
		int i = *a;
		int j = *b;
		if (abs(B[i][j]) > eps) {
			T = (B[j][j] - B[i][i]) / (2. * B[i][j]);
			t = sgn(T) / (abs(T) + sqrt(1. + T * T));
			c = 1. / sqrt(1. + t * t);
			s = c * t;
		}
		else {
			c = 1.;
			s = 0;
		}
		A_cur[i][i] = c * c * B[i][i] - 2 * s * c * B[i][j] + s * s * B[j][j];
		A_cur[j][j] = s * s * B[i][i] + 2 * s * c * B[i][j] + c * c * B[j][j];
		A_cur[i][j] = (c * c - s * s) * B[i][j] + s * c * (B[i][i] - B[j][j]);
		A_cur[j][i] = A_cur[i][j];
		for (int k = 0; k < n; ++k) {
			if ((k != i) && (k != j)) {
				A_cur[i][k] = c * B[i][k] - s * B[j][k];
				A_cur[k][i] = A_cur[i][k];
				A_cur[j][k] = s * B[i][k] + c * B[j][k];
				A_cur[k][j] = A_cur[j][k];
			}
		}
		for (int k = 0; k < n; ++k) {
			for (int l = 0; l < n; ++l) {
				if ((k != i) && (k != j) && (l != i) && (l != j)) {
					A_cur[k][l] = B[k][l];
				}
			}
		}
		error = check_conv(A_cur);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				B[i][j] = A_cur[i][j];
			}
		}
	} while (error >= eps);
	print_matrix(B);
	cout << "Error: " << setw(15) << right << scientific << error << endl;
	delete a;
	delete b;
	return B[0][0];
}

double sc_eigenv_max(double(*A)[n]) { //max by abs eigenvalue by scalar product method
	double max = 0;
	double x[n] = { 0 };
	double y[n] = { 0 };
	for (int i = 0; i < n; ++i) {
		y[i] = 1;
	}
	double s, norm, s_cur, t, norm_cur, max_cur, error;
	s = prod(y, y); //null iter
	norm = sqrt(s);
	for (int i = 0; i < n; ++i) {
		y[i] /= norm; //normalized x0
	}
	do {
		double y_cur[n] = { 0 };
		//first iter
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				y_cur[i] += A[i][j] * y[j];
			}
		}
		s_cur = prod(y_cur, y_cur);
		t = prod(y_cur, y);
		norm_cur = sqrt(s_cur);
		max_cur = s_cur / t;
		error = abs(max_cur - max);
		max = max_cur;
		for (int i = 0; i < n; ++i) {
			y[i] = y_cur[i] / norm_cur; //normalized x0
		}
	} while (error >= eps);
	return max;
}

void matrix_product(double(*Res)[n], double(*A)[n], double(*B)[n]) {
	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Res[i][k] += A[i][j] * B[j][k];
			}
		}
	}
}

double min_by_abs_eigenvalue(double(*A)[n]) {
	double C[n][n] = { 0 };
	double max = sc_eigenv_max(A);
	matrix_product(C, A, A);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			C[i][j] *= -1. / (max * max);
		}
		C[i][i] += 1;
	}

	double max_C = sc_eigenv_max(C);
	return  sqrt(1 - max_C) * max;

}

double min_eigenvalue(double(*A)[n]) {
	double B[n][n] = { 0 };
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B[i][j] -= A[i][j];
		}
		B[i][i] += sc_eigenv_max(A);
	}
	return sc_eigenv_max(A) - sc_eigenv_max(B);
}



int main() {
	std::cout.precision(8);
	//double A[n][n] = { {4,-30,60,-35},{-30,300,-675,420},{60,-675,1620,-1050},{-35,420,-1050,700} };
	double A[n][n] = { 0 };
	for (int i = 1; i < n + 1; i++) {
		A[i - 1][i - 1] = 1. / i + i * n * n;
		for (int j = i + 1; j < n + 1; j++) {
			A[i - 1][j - 1] = pow(-1, j) * j;
			A[j - 1][i - 1] = A[i - 1][j - 1];
		}
	}

	cout << "Input matrix A:" << endl;
	print_matrix(A);
	cout << "SP alg for A" << endl;
	cout << setw(15) << right << scientific << "Max by absolute value eigenvalue: " << sc_eigenv_max(A) << endl;
	cout << setw(15) << right << scientific << "Min by absolute value eigenvalue: " << min_by_abs_eigenvalue(A) << endl;
	cout << setw(15) << right << scientific << "Min eigenvalue: " << min_eigenvalue(A) << endl;

	cout << "Jacoby for A" << endl;
	Jacoby_eigenv(A);
	//cout << "Jacoby for C" << endl;
	//Jacoby_eigenv(C);

	system("pause");
}