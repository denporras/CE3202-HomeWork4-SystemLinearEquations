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

//+++Test cases+++
//testLU
TEST(testLU, NormError) {

	Matrix<double> LU;
	Matrix<double> M = {{4,5,7},{1,2,4},{1,1,3}};

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();
	TestDescomposition<double> * test = new TestDescomposition<double>();

	d->lu(M,LU);
	double norm = test->testLU(M,LU);

    ASSERT_EQ(norm , 0);

}

//testQR
TEST(testQR, NormError) {

	Matrix<double> Q;
	Matrix<double> R;
	Matrix<double> M = {{4,5,7},{1,2,4},{1,1,3}};

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();
	TestDescomposition<double> * test = new TestDescomposition<double>();

	d->qr(M,Q,R);
	double norm = test->testQR(M,Q,R);

    ASSERT_EQ(norm , 0);

}

//testInvert (Must be passed)
TEST(testInvert, NoThrow) {

	Matrix<double> M = {{4,5,7},{1,2,4},{1,1,3}};
	Matrix<double> Mi;

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();

    ASSERT_NO_THROW(d->inverse(M,Mi));

}

//testInvert (Must fail)
TEST(testInvert, Throw) {

	Matrix<double> M = {{0,0,0},   //Normal 1
						{0,0,0},
						{0,0,0}};
	Matrix<double> Mi;

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();

	ASSERT_ANY_THROW(d->inverse(M,Mi));


}

//testMultiplication (Must be passed)
TEST(Multiplication, NoThrow) {

	Matrix<double> M = {{4,5,7},{1,2,4},{1,1,3}}; ;


    ASSERT_NO_THROW(M*M);

}

//testMultiplication (Must fail)
TEST(Multiplication, Throw) {

	Matrix<double> M = {{4,5,7},{1,2,4},{1,1,3}}; ;

	Matrix<double> M2 = {{1,1},
						{1,2}};

	ASSERT_ANY_THROW(M*M2);

}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
