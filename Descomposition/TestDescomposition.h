/**
 * @file TestDescomposition.h
 * @brief Class that implements tests to LU and QR methods
 * @author Kevin Alfaro
 * @date 27 de sept. de 2017
 */

#ifndef DESCOMPOSITION_TESTDESCOMPOSITION_H_
#define DESCOMPOSITION_TESTDESCOMPOSITION_H_

#include "../Matrix/Matrix.hpp"
#include <cmath>
#include <gtest/gtest.h>

template<typename T>
class TestDescomposition {
public:
	TestDescomposition();
	T testLU(anpi::Matrix<T> &A, anpi::Matrix<T> &LU);
	T testQR(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q,anpi::Matrix<T> &R);

private:
	void getUpperLowerTringular(anpi::Matrix<T> &M,anpi::Matrix<T> &UT,anpi::Matrix<T> &LT);
	T getNorm(anpi::Matrix<T> &M);
};

/**
 * @brief Constructor by default.
 */
template<typename T>
TestDescomposition<T>::TestDescomposition() {}

/**
 * @brief Test LU method comparing original A and A-from-LU reconstruction
 * @param A: Original matrix
 * @param LU: Matrix from LU method
 * @return Norm  of Difference between  A and reconstructed A
 */
template<typename T>
T TestDescomposition<T>::testLU(anpi::Matrix<T> &A, anpi::Matrix<T> &LU){

	ASSERT_EQ(A.rows(),A.cols());   //A matrix is square
	ASSERT_EQ(LU.rows(),LU.cols()); //LU matrix is square

	anpi::Matrix<T> LT(LU.rows(),LU.cols(),T(0));
	anpi::Matrix<T> UT(LU.rows(),LU.cols(),T(0));
	anpi::Matrix<T> A_aux(A.rows(),A.cols(),T(0));

	//Separate LU to L and U
	this->getUpperLowerTringular(LU,UT,LT);

	A_aux = LT*UT; //Reconstruct A from LU
	anpi::Matrix<T> A_difference = A-A_aux; //Difference between original A and reconstructed-from-LU A

	T norm = this->getNorm(A_difference);
	ASSERT_LT(norm,0.01);//Error less than 0.01%


	return norm;

}

/**
 * @brief Test QR method comparing original A and A-from-LU reconstruction
 * @param A: Original matrix
 * @param Q: Matrix Q from QR method
 * @param R: Matrix R from QR method
 * @return Norm  of Difference between  A and reconstructed A
 */
template<typename T>
T TestDescomposition<T>::testQR(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q, anpi::Matrix<T> &R){

	ASSERT_EQ(A.rows(),A.cols()); //A matrix is square
	ASSERT_EQ(Q.rows(),Q.cols()); //Q matrix is square
	ASSERT_EQ(R.rows(),R.cols()); //R matrix is square

	anpi::Matrix<T> A_aux(A.rows(),A.cols(),T(0));


	A_aux = Q*R; //Reconstruct A from QR
	anpi::Matrix<T> A_difference = A-A_aux; //Difference between original A and reconstructed-from-QR A

	T norm = this->getNorm(A_difference);
	ASSERT_LT(norm,0.01);//Error less than 0.01%

	return norm;

}

/**
 * @brief Substract Upper and Lower triangular matrixs
 * @param M: Original matrix
 * @param LT: Destination Lower triangular 0 matrix with the same dimensions than M
 * @param UT: Destination Upper triangular 0 matrix with the same dimensions than M
 */
template<typename T>
void TestDescomposition<T>::getUpperLowerTringular(anpi::Matrix<T> &M,anpi::Matrix<T> &UT,anpi::Matrix<T> &LT){

	//Test if original matrix and destination matrix have the same size and they are square matrix
	ASSERT_EQ(M.rows(),M.cols());  //Square
	ASSERT_EQ(M.rows(),LT.rows()); //M rows = LT rows
	ASSERT_EQ(M.cols(),LT.cols()); //M rows = LT cols
	ASSERT_EQ(M.rows(),UT.rows()); //M rows = UT rows
	ASSERT_EQ(M.cols(),UT.cols()); //M rows = UT cols

	//Upper triangular matrix
	for(int i = 0; i < M.rows(); i++){
		for(int j = i ; j < M.cols(); j++){
			UT[i][j] = M[i][j];
		}
	}

	//Lower Triangular matrix
	for(int i = 0; i < M.rows(); i++){
		for(int j = 0; j <= i; j++){
			if(i == j){
				LT[i][j] = 1; //Diagonal values
			}
			else{
				LT[i][j] = M[i][j];
			}
		}
	}

}

/**
 * @brief Calculate the matrix norm
 * @param M: Matrix
 * @return Norm  of M
 */
template<typename T>
T TestDescomposition<T>::getNorm(anpi::Matrix<T> &M){

	T a_ij = T(0);
	T norm = T(0);

	for(int i = 0; i < M.rows(); i++){
		for(int j = 0; j < M.rows(); j++){
			a_ij+=abs(M[i][j]*M[i][j]);
		}
	}

	norm = sqrt(a_ij);

	return norm;
}

#endif /* DESCOMPOSITION_TESTDESCOMPOSITION_H_ */
