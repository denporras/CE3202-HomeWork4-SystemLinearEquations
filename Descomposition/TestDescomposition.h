/*
 * TestDescomposition.h
 *
 *  Created on: 27 de set. de 2017
 *      Author: kevin
 */

#ifndef DESCOMPOSITION_TESTDESCOMPOSITION_H_
#define DESCOMPOSITION_TESTDESCOMPOSITION_H_

//#include <boost/test/unit_test.hpp>
#include "../Matrix/Matrix.hpp"

template<typename T>
class TestDescomposition {
public:
	TestDescomposition();
	T testLU(anpi::Matrix<T> &A, anpi::Matrix<T> &LU);
	T testQR(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q,anpi::Matrix<T> &R);

private:
	void getUpperLowerTringular(anpi::Matrix<T> &M,anpi::Matrix<T> &UT,anpi::Matrix<T> &LT);
};

template<typename T>
TestDescomposition<T>::TestDescomposition() {

}

template<typename T>
T TestDescomposition<T>::testLU(anpi::Matrix<T> &A, anpi::Matrix<T> &LU){
	anpi::Matrix<T> LT(LU.rows(),LU.cols(),T(0));
	anpi::Matrix<T> UT(LU.rows(),LU.cols(),T(0));
	this->getUpperLowerTringular(LU,UT,LT);
	A = LT*UT;
}

template<typename T>
T TestDescomposition<T>::testQR(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q, anpi::Matrix<T> &R){

}

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


#endif /* DESCOMPOSITION_TESTDESCOMPOSITION_H_ */
