/*#include <iostream>
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

//Test cases
TEST(testLU, NormError) {

	Matrix<double> LU;
	Matrix<double> M = {{1,1,3},   //Normal 1
						 {1,2,4},
						 {4,5,7}};

	MatrixDescomposition<double> * d = new MatrixDescomposition<double>();
	TestDescomposition<double> * test = new TestDescomposition<double>();

	d->lu(M,LU);
	double norm = test->testLU(M,LU);

    ASSERT_LT(norm , 0.01);

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}*/
