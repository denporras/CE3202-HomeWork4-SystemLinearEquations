//============================================================================
// Name        : Resolucion de Sistemas de Ecuacione Lineales
// Author      : Kevin Alfaro Vega, Dennis Porras Barrantes, David Gómez Vargas
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

using namespace std;
using namespace anpi;

template <typename T>
void printMatrix(Matrix<T> &m){
	for(int i = 0; i < m.rows(); i++){
		cout << "|";
		for(int j = 0; j < m.cols(); j++)
			cout << "[" << m[i][j] << "]\t";
		cout << "|" <<endl;
	}
}


int main() {
	Matrix<double> A= {{25,5,1},{64,8,1},{144,12,1}};
	Matrix<double> LU;
	LUDescomposition<double> * d = new LUDescomposition<double>();
	d->lu(A, LU);
	printMatrix(LU);
	return 0;
}
