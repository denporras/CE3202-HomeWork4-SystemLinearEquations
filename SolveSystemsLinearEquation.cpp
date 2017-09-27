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
#include "Descomposition/TestDescomposition.h"
#include "Matrix/Matrix.hpp"

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

int main()
{
	Matrix<double> A= {{25,5,1},{64,8,1},{144,12,1}};
	Matrix<double> B= {{25,5,1},{64,8,1},{144,12,1}};

	Matrix<double> C= {{1,2,3}};
	Matrix<double> D= {{1,2},{1,2},{4,4}};

	Matrix<double> RES(C.rows(),D.cols(),1);


	Matrix<double> M= {{25,5,1},{64,8,1},{144,12,1}};
	Matrix<double> LT(A.rows(),A.cols(),double(0));
	Matrix<double> UT(A.rows(),A.cols(),double(0));

	TestDescomposition<double> * test = new TestDescomposition<double>();
	test->getUpperLowerTringular(M,UT,LT);

	cout << "Upper";
	printMatrix(UT);
	cout << "Lower";
	printMatrix(LT);

	/*try{
		RES = C*D;
		printMatrix(RES);
	}
	catch (const std::exception& error){
		std::cerr << "Exception: " << error.what() << std::endl;
	}
	catch (...){

		std::cerr << "Exception: unknown" << std::endl;
	}*/

	return 0;
}

/*int main() {
	Matrix<double> A= {{25,5,1},{64,8,1},{144,12,1}};
	Matrix<double> LU;
	LUDescomposition<double> * d = new LUDescomposition<double>();
	d->lu(A, LU);
	printMatrix(LU);
	return 0;
}*/
