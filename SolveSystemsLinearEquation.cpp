//============================================================================
// Name        : Resolucion de Sistemas de Ecuacione Lineales
// Author      : Kelvin Alfaro Vega, Dennis Porras Barrantes, David Gómez Vargas
// Version     : 1.0
// Copyright   : Assignment for the course Numerical Analysis of the Costa Rica Institute of Technology
// Description : Program that solves systems of linear equations.
//============================================================================

#include <iostream>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstdlib>

#include "Descomposition/LUDescomposition.h"
#include "Descomposition/QRDescomposition.h"

using namespace std;
using namespace anpi;

template <typename T>
void printMatrix(anpi::Matrix<T> &m){
	for(int i = 0; i < m.rows(); i++){
		cout << "|\t";
		for(int j = 0; j < m.cols(); j++)
			cout << "[" << m[i][j] << "]\t";
		cout << "|" <<endl;
	}
}

template <typename T>
void printVector(vector<T> &v){
	cout << "[\t";
	for(int i = 0; i < v.size(); i++){
		cout << v.at(i) << "\t";
	}
	cout << "]" << endl;
}

int main() {
	Matrix<double> A= {{1,1,3},{1,2,4},{4,5,7}};
	LUDescomposition<double> * d = new LUDescomposition<double>();
	vector<double> x;
	vector<double> b = {{7,11,20}};
	d->solve(A, x, b);
	printVector(x);

	Matrix<double> A_1= {{2, 1, 3, 5},{-1, 0, 7, 1},{0, -1, -1, 3},{-3, 7, 4, 3},{1, 6, 4, -3}};
	Matrix<double> Q;
	Matrix<double> R;
	QRDescomposition<double> * des = new QRDescomposition<double>();
	des->qr(A_1,Q,R);
	printMatrix(Q);
	printMatrix(R);
	return 0;
}
