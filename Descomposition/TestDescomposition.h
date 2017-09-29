/**
 * @file TestDescomposition.h
 * @brief Class that implements tests to LU and QR methods
 * @author Kevin Alfaro
 * @date 27 de sept. de 2017
 */

#ifndef DESCOMPOSITION_TESTDESCOMPOSITION_H_
#define DESCOMPOSITION_TESTDESCOMPOSITION_H_


#include <cmath>
#include "../Matrix/Matrix.hpp"

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


	anpi::Matrix<T> LT(LU.rows(),LU.cols(),T(0));
	anpi::Matrix<T> UT(LU.rows(),LU.cols(),T(0));
	anpi::Matrix<T> A_aux(A.rows(),A.cols(),T(0));

	//Separate LU to L and U
	this->getUpperLowerTringular(LU,UT,LT);

	A_aux = LT*UT; //Reconstruct A from LU
	anpi::Matrix<T> A_difference = A-A_aux; //Difference between original A and reconstructed-from-LU A

	T norm = this->getNorm(A_difference);

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


	anpi::Matrix<T> A_aux(A.rows(),A.cols(),T(0));


	A_aux = Q*R; //Reconstruct A from QR
	anpi::Matrix<T> A_difference = A-A_aux; //Difference between original A and reconstructed-from-QR A

	T norm = this->getNorm(A_difference);

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
