//============================================================================
// Name        : Resolucion de Sistemas de Ecuacione Lineales
// Author      : Kevin Alfaro Vega, Dennia Porras Barrantes, David Gómez Vargas
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

#include "Descomposition/TestDescomposition.h"
#include "Descomposition/MatrixDescomposition.h"
#include <gtest/gtest.h>
#include "Matrix/Matrix.hpp"

using namespace std;
using namespace anpi;


/**
 * @brief Print a matrix
 * @param m: Matrix to print
 */
template <typename T>
void printMatrix(anpi::Matrix<T> &m){
	for(int i = 0; i < m.rows(); i++){
		cout << "|\t";
		for(int j = 0; j < m.cols(); j++)
			cout << "[" << m[i][j] << "]\t";
		cout << "|" <<endl;
	}
}

/**
 * @brief Print a vector
 * @param v: Vector to print
 */
template <typename T>
void printVector(vector<T> &v){
	cout << "[\t";
	for(int i = 0; i < v.size(); i++){
		cout << v.at(i) << "\t";
	}
	cout << "]" << endl;
}

template<class T>
void setUp(int m, int f, int s) {

	//Matrix
	Matrix<T> M_s;//Selected M
	Matrix<T> Q;
	Matrix<T> R;
	Matrix<T> Mi;
	Matrix<T> LU;

	Matrix<T> M_1 = {{1,1,3},   //Normal 1
					 {1,2,4},
					 {4,5,7}};

	Matrix<T> M_2 = {{8,14,3},  //Normal 2
					 {3,2,7},
					 {6,1,9}};

	Matrix<T> M_3= {{1,0.01},   //Mal cond 1
					{0.99,1}};

	Matrix<T> M_4= {{4,  5 },   //Mal cond 2
					{4.1,5}};

	Matrix<T> M_5= {{1,1,1},    //Non invertible
					{1,1,1},
					{1,1,1}};


	//Equation systems
	vector<T> x; //Result vector
	Matrix<T> A_s;//Selected A
	vector<T> b_s;//Selected b

	Matrix<T> A_1= {{1,1,3},
					{1,2,4},
					{4,5,7}};

	vector<T> b_1 = {{7,11,20}};

	Matrix<T> A_2 = {{1,0,5,7,0,7},
					{3,3,44,6,8,8},
					{6,4,3,3,4,6 },
					{89,8,7,6,5,4},
					{5,6,8,8,9,8 },
					{65,4,3,3,4,5}};

	vector<T> b_2 = {{4,7,7,4,7,7}};

	Matrix<T> A_3= {{2,1,3,5  }, //Wrong dimensions
					{-1,0,7,1 },
					{0,-1,-1,3},
					{-3,7,4,3 },
					{1,6,4,-3 }};

	vector<T> b_3 = {{7,11,20}};
	T norm;

	//Select equation system
	switch(s){
	case 1:
		A_s = A_1;
		b_s = b_1;

		break;
	case 2:
		A_s = A_2;
		b_s = b_2;
		break;
	case 3:
		A_s = A_3;
		b_s = b_3;
		break;
	default:
		A_s = A_1;
		b_s = b_1;

	}

	//Select matrix
	switch(m){
	case 1:
		M_s = M_1;

		break;
	case 2:
		M_s = M_2;
		break;
	case 3:
		M_s = M_3;
		break;
	case 4:
		M_s = M_4;
		break;
	case 5:
		M_s = M_5;
		break;
	default:
		M_s = M_1;

	}

	MatrixDescomposition<T> * d = new MatrixDescomposition<T>();
	TestDescomposition<T> * test = new TestDescomposition<T>();

	//Select function
	switch(f){
	case 1://lu

		d->lu(M_s,LU);
		cout << "LU Matrix" << endl;
		printMatrix(LU);

		break;
	case 2: //qr
		d->qr(M_s,Q,R);
		cout << "Q Matrix" << endl;
		printMatrix(Q);
		cout << "R Matrix" << endl;
		printMatrix(R);

		break;
	case 3: //testLU

		d->lu(M_s,LU);
		norm = test->testLU(M_s,LU);

		cout << "Norm of difference between original matrix and reconstructed-from-LU matrix:"<<endl;
		cout << norm <<endl;

		break;
	case 4: //testQR

		d->lu(M_s,LU);
		norm = test->testQR(M_s,Q,R);

		cout << "Norm of difference between original matrix and reconstructed-from-QR matrix:"<<endl;
		cout << norm <<endl;


		break;
	case 5: //solveLU

		d->solveLU(A_s,x,b_s);
		cout << "Vector x (solution)" << endl;
		printVector(x);


		break;
	case 6: //solveQR

		d->solveQR(A_s,x,b_s);
		cout << "Vector x (solution)" << endl;
		printVector(x);

			break;
	case 7: //invert
		d->inverse(M_s,Mi);
		cout << "Inverted matrix" << endl;
		printMatrix(Mi);

			break;
	default:
		d->lu(M_s,LU);
		cout << "LU Matrix" << endl;
		printMatrix(LU);
	}
}


int main() {


	int f, p, m,s;


	cout << "Choose a number (1,2) to select the precision:" << endl
			<< "1. Float" << endl << "2. Double" << endl;
	cin >> p;

	cout << "Choose a number (1,2,3,4,5,6,7) to select the function:" << endl;
	cout << "1. lu" << endl;
	cout << "2. qr" << endl;
	cout << "3. testLU" << endl;
	cout << "4. testQR" << endl;
	cout << "5. solveLU" << endl;
	cout << "6. solveQR" << endl;
	cout << "7. invert" << endl;
	cin >> f;

	if(f == 5 or f == 6){
		cout << "Choose a number (1,2,3) to select the option Ax = b system: (see README)" << endl;
		cin >> s;

	}
	else{
		cout << "Choose a number (1,2,3,4,5) to select the matrix: (see README)" << endl;
		cin >> m;
	}

	if (p == 1) {
		setUp<float>(m, f, s);
	} else if (p == 2) {
		setUp<double>(m, f ,s);
	}

	return 0;
}



