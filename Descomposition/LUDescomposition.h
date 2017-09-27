/**
 * @brief Class that applies the LU descomposition
 * @author Dennis Porras Barrantes
 * @date 23/09/2017
 */

#ifndef DESCOMPOSITION_LUDESCOMPOSITION_H_
#define DESCOMPOSITION_LUDESCOMPOSITION_H_

#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include "../Matrix/Matrix.hpp"

using namespace std;

template<typename T>
class LUDescomposition {
public:
	LUDescomposition();
	void lu(const anpi::Matrix<T> &A, anpi::Matrix<T> &LU);
	void solve(const anpi::Matrix<T> &A, vector<T> &x, const vector<T> &b);
	void inverse(const anpi::Matrix<T> &A ,anpi::Matrix<T> &Ai);
private:
	int n;
	vector<T> index;
	T determinant;
	anpi::Matrix<double> luMatrix;
};

template<typename T>
LUDescomposition<T>::LUDescomposition() {
	this->n = 0;
	this->determinant = T(0);
}

template<typename T>
void LUDescomposition<T>::lu(const anpi::Matrix<T> &A, anpi::Matrix<T> &LU){
	this->n = A.rows();
	LU = A;

	const T SMALL = 1.0e-40;
	int i, i_max, j, k;
	T big, tmp;
	vector<T> scaling;
	this->determinant = T(0);
	for(i = 0; i < this->n; i++){
		big = T(0);
		for(j = 0; j < this->n; j++){
			if((tmp = abs(LU[i][j])) > big)
				big = tmp;
		}
		if(abs(big) < numeric_limits<T>::epsilon()){
			throw runtime_error("Singular matrix in method: \n void lu(const Matrix<T>& A, Matrix<T>& LU)");
		}
		scaling.push_back(T(1)/big);
	}

	for(k = 0; k < this->n; k++){
		big = T(0);
		for(i = k; i < this->n; i++){
			tmp = scaling.at(i) * abs(LU[i][k]);
			if(tmp > big){
				big = tmp;
				i_max = i;
			}
		}
		if(k != i_max){
			for(j = 0; j < this->n; j++){
				tmp = LU[i_max][j];
				LU[i_max][j] = LU[k][j];
				LU[k][j] = tmp;
			}
			this->determinant *= -1;
			scaling.at(i_max) = scaling.at(k);
		}
		this->index.push_back(i_max);
		if(abs(LU[k][k]) < numeric_limits<T>::epsilon())
			LU[k][k] = SMALL;
		for(i = k+1;i < this->n; i++){
			tmp = LU[i][k] /= LU[k][k];
			for(j= k+1; j < this->n; j++)
				LU[i][j] -= tmp*LU[k][j];
		}

	}
}

template<typename T>
void LUDescomposition<T>::solve(const anpi::Matrix<T> &A, vector<T> &x, const vector<T> &b) {
	int i, ip, j;
	int ii = 0;
	T sum;
	this->lu(A, this->luMatrix);
	if(b.size() != this->n)
		throw runtime_error("The rows dimension is not correct on void LUDescomposition<T>::solve(vector<T>& b, vector<T>& x)");
	x = b;
	for(i = 0; i < this->n; i++){
		ip = this->index.at(i);
		sum = x.at(ip);
		x.at(ip) = x.at(i);
		if(ii != 0){
			for(j = ii-1; j < i; j++){
				sum -= (this->luMatrix[i][j])*(x.at(j));
			}
		}else if(abs(sum) > numeric_limits<T>::epsilon())
			ii = i+1;
		x.at(i) = sum;
	}for(i = n-1; i >= 0; i--){
		sum = x.at(i);
		for(j = i+1; j < n; j++)
			sum -= this->luMatrix[i][j]*x.at(j);
		x.at(i) = sum/(this->luMatrix[i][i]);
	}
}

template<typename T>
void LUDescomposition<T>::inverse(const anpi::Matrix<T> &A, anpi::Matrix<T>& Ai) {

}

#endif /* DESCOMPOSITION_LUDESCOMPOSITION_H_ */
