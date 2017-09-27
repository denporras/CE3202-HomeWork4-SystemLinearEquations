/*
 * MatrixDescomposition.h
 *
 *  Created on: 27 sept. 2017
 *      Author: dennis
 */

#ifndef DESCOMPOSITION_MATRIXDESCOMPOSITION_H_
#define DESCOMPOSITION_MATRIXDESCOMPOSITION_H_

#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include "../Matrix/Matrix.hpp"

using namespace std;

namespace anpi {

template<typename T>
class MatrixDescomposition {
public:
	MatrixDescomposition();
	void lu(const Matrix<T> &A, Matrix<T> &LU);
	bool solveLU(const Matrix<T> &A, vector<T> &x, const vector<T> &b);
	void inverse(const Matrix<T> &A ,Matrix<T> &Ai);
	void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
    bool solveQR(const Matrix<T> &A, vector<T> &x, const vector<T> &b);
private:
	Matrix<T> ext_prod(Matrix<T> &a, Matrix<T> &b);
    Matrix<T> scal_mat(Matrix<T> &a, int b);
    Matrix<T> transpose(Matrix<T> &a);
	void inverse_aux(const Matrix<T> &B, Matrix<T> &X);
	bool inverseFlag;
	int n;
	vector<T> index;
	Matrix<double> luMatrix;
};



template<typename T>
inline MatrixDescomposition<T>::MatrixDescomposition() {
	this->n = 0;
	this->inverseFlag = false;
}

template<typename T>
inline void MatrixDescomposition<T>::lu(const Matrix<T>& A,
		Matrix<T>& LU) {
	this->n = A.rows();
	if(this->n != A.cols())
		throw runtime_error("'A' matrix is not square in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
	this->index.clear();
	LU = A;

	const T SMALL = 1.0e-40;
	int i, i_max, j, k;
	T big, tmp;
	vector<T> scaling;
	for(i = 0; i < this->n; i++){
		big = T(0);
		for(j = 0; j < this->n; j++){
			if((tmp = abs(LU[i][j])) > big)
				big = tmp;
		}
		if(abs(big) < numeric_limits<T>::epsilon()){
			throw runtime_error("Singular matrix in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
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
inline bool MatrixDescomposition<T>::solveLU(const Matrix<T>& A,
		vector<T>& x, const vector<T>& b) {
	bool result = true;
	int i, ip, j;
	int ii = 0;
	T sum;
	if(!this->inverseFlag)
		this->lu(A, this->luMatrix);

	if(b.size() != this->n)
		throw runtime_error("The rows dimension is not correct in void LUDescomposition<T>::solve(vector<T>& b, vector<T>& x)");

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
	return result;
}

template<typename T>
inline void MatrixDescomposition<T>::inverse(const Matrix<T>& A,
		Matrix<T>& Ai) {
	this->lu(A,this->luMatrix);
	this->inverseFlag = true;
	Ai = anpi::Matrix<T>(this->n, this->n, T(0));
	for(int i = 0; i < this->n; i++){
		Ai[i][i] = 1;
	}
	this->inverse_aux(Ai,Ai);
	this->inverseFlag = false;
}

template<typename T>
inline void MatrixDescomposition<T>::inverse_aux(const Matrix<T>& B,
		Matrix<T>& X) {
	int i, j;
	int m = B.cols();
	if(B.rows() != this->n || X.rows() != this->n || B.cols() != X.cols()){
		throw runtime_error("Bad sizes for the matrix in void LUDescomposition<T>::inverse_aux(const anpi::Matrix<T>& B, anpi::Matrix<T>& X)");
	}
	vector<T> tmp;
	for(j = 0; j < m; j++){
		tmp.clear();
		for(i = 0; i < this->n; i++)
			tmp.push_back(B[i][j]);
		this->solveLU(B, tmp, tmp);
		for(i = 0; i < this->n; i++)
			X[i][j] = tmp.at(i);
	}
}


template<typename T>
void MatrixDescomposition<T>::qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R){
	int m = A.rows();
	int n = A.cols();
    T mag, alpha;
    T x = 0;
	Matrix<T> u(m,1,x);
	Matrix<T> v(m,1,x);
	Matrix<T> P(m,m,x);
    for (int i = 0; i<m;++i){
        P(i,i) = 1;
    }
	Matrix<T> I = P;
    Q = P;
    R = A;
    for (int i = 0; i < n; ++i) {         
        u.fill(x);
        v.fill(x);
        mag = 0.0;
        for (int j = i; j < m; ++j) {
            u(j,0) = R(j,i);
            mag += u(j,0) * u(j,0);
        }
        mag = sqrt(mag);
        alpha = u(i,0) < 0 ? mag : -mag;
        mag = 0.0;
        for (int j = i; j < m; j++) {
            v(j,0) = j == i ? u(j,0) + alpha : u(j,0);
            mag += v(j,0) * v(j,0);
        }
        mag = sqrt(mag);
        if (mag < 0.0000000001) continue;
        for (int j = i; j < m; j++) v(j,0) /= mag;   
        Matrix<T> t(1,m,x);
        for (int j = 0; j < m; ++j){
            t(0,j) = v(j,0);
        }
        Matrix<T> w = ext_prod(v,t);
        P = I - scal_mat(w,2);
        R = P*R;
        Q = Q*P;
    }
}


template <typename T>
Matrix<T> MatrixDescomposition<T>::ext_prod(Matrix<T> &a, Matrix<T> &b){
     Matrix<T> r(a.rows(),b.cols(),1);
    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < b.cols(); j++)
            r(i,j) = a(i,0)*b(0,j);
    }
    return r;
}


template <typename T>
Matrix<T> MatrixDescomposition<T>::scal_mat(Matrix<T> &a, int b){
    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            a(i,j) = a(i,j)*b;
    }
    return a;
}

template <typename T>
Matrix<T> MatrixDescomposition<T>::transpose(Matrix<T> &a){
    T x = 0;
    Matrix<T> r(a.cols(),a.rows(),x);
    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            r(j,i) = a(i,j);
    }
    return r;
}


template <typename T>
bool MatrixDescomposition<T>::solveQR(const Matrix<T> &A, vector<T> &x, const vector<T> &b){
    vector<T> bp;
    Matrix<double> Q;
    Matrix<double> R;
    qr(A,Q,R);
    int n = b.size();
    for (int i =0; i< n;++i){
        x.push_back(0);
    }

    Matrix<double> Qt = transpose(Q);
    for (int i = 0; i< Qt.rows();++i){
       bp.push_back(0);
       for (int j = 0; j<Qt.cols();++j){
        bp[i] += Qt(i,j) * b[j];
        }
    }

    x[n-1] = bp[n-1] / R[n-1][n-1];
    for (int i = n-2; i >= 0; --i){
        T sum=bp[i];
        for (int j = i+1;j < n; ++j){
            sum -= R[i][j] * x[j];
        }
        x[i]=sum / R[i][i];
    }
    return 0;

}

} /* namespace anpi */
#endif /* DESCOMPOSITION_MATRIXDESCOMPOSITION_H_ */