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

#include "Descomposition/TestDescomposition.h"
#include "Matrix/Matrix.hpp"
#include "Descomposition/MatrixDescomposition.h"


using namespace std;
using namespace anpi;



template<class T>
void setUp(int m, int f, int s) {

	//Matrix
	Matrix<T> M_s;//Selected M
	Matrix<T> LU;
	Matrix<T> Q;
	Matrix<T> R;
	Matrix<T> Mi;

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
	vector<T> A_s;//Selected A
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


}


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

void setUp(int m, int f, int s) {

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




/**
	Matrix<double> A= {{1,1,3},{1,2,4},{4,5,7}};
	vector<double> x;
	vector<double> b = {{7,11,20}};
	LUDescomposition<double> * d = new LUDescomposition<double>();
	d->solve(A, x, b);
	printVector(x);

	vector<double> x_1;
	vector<double> b_1 = {{7,11,20}};
	Matrix<double> A_1= {{2, 1, 3, 5},{-1, 0, 7, 1},{0, -1, -1, 3},{-3, 7, 4, 3},{1, 6, 4, -3}};
	QRDescomposition<double> * des = new QRDescomposition<double>();
	//des->qr(A_1,Q,R);
	des->solveQR(A,x_1,b_1);
	printVector(x_1);
*/


/*	Matrix<double> A = {{1,0,5,7,0,7},
						{3,3,44,6,8,8},
						{6,4,3,3,4,6},
						{89,8,7,6,5,4},
						{5,6,8,8,9,8},
						{65,4,3,3,4,5}};
	vector<double> x;
	vector<double> b = {{4,7,7,4,7,7}};
	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();
	d->solveLU(A, x, b);
	printVector(x);

	cout << endl;

	Matrix<double> M = {{5,6,4},
						{5,8,8},
						{5,1,3}};
	Matrix<double> Mi;
	d->inverse(M, Mi);
	//printMatrix(Mi);
	Matrix<double> A_1= {{1,0,5,7,0,7},
						{3,3,44,6,8,8},
						{6,4,3,3,4,6},
						{89,8,7,6,5,4},
						{5,6,8,8,9,8},
						{65,4,3,3,4,5}};
	vector<double> x_1;
	vector<double> b_1 = {{4,7,7,4,7,7}};
	d->solveQR(A_1,x_1,b_1);
	printVector(x_1);*/

	/*

	TestDescomposition<double> * test = new TestDescomposition<double>();

	Matrix<double> M= {{1,1,3},{1,2,4},{4,5,7}};

	Matrix<double> LU= {{0,0,0},{0,0,0},{0,0,0}};

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();
	d->lu(M,LU);

	cout << "LU"<<endl;
	printMatrix(LU);

	double err = test->testLU(M,LU);

	cout << "M"<<endl;
	printMatrix(M);

	cout << "err"<<endl;
	cout << err<<endl;*/


	/*try{
			a = b*c;
		}
		catch (const std::exception& error){
			std::cerr << "Exception: " << error.what() << std::endl;
		}*/




	return 0;
}


