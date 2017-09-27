/**
 * @brief Class that applies the LU descomposition
 * @author Dennis Porras Barrantes
 * @date 23/09/2017
 */

#ifndef DESCOMPOSITION_QRDESCOMPOSITION_H_
#define DESCOMPOSITION_QRDESCOMPOSITION_H_

#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include "../Matrix/Matrix.hpp"

using namespace std;
using namespace anpi;

template<typename T>
class QRDescomposition {
public:
	QRDescomposition();
	void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
private:
    Matrix<T> ext_prod(Matrix<T> &a, Matrix<T> &b);
    Matrix<T> mat_prod(Matrix<T> &a, Matrix<T> &b);
    Matrix<T> scal_mat(Matrix<T> &a, int b);
};

template<typename T>
QRDescomposition<T>::QRDescomposition() {
}

template<typename T>
void QRDescomposition<T>::qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R){
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
        R = mat_prod(P,R);
        Q = mat_prod(Q,P);
    }
}



template <typename T>
Matrix<T> QRDescomposition<T>::ext_prod(Matrix<T> &a, Matrix<T> &b){
     Matrix<T> r(a.rows(),b.cols(),1);
    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < b.cols(); j++)
            r(i,j) = a(i,0)*b(0,j);
    }
    return r;
}

template <typename T>
Matrix<T> QRDescomposition<T>::mat_prod(Matrix<T> &a, Matrix<T> &b){
     T x = 0;
     Matrix<T> r(a.rows(),b.cols(),x);
     for(int i=0; i<a.rows(); i++)
        for(int j=0; j<b.cols(); j++)
            for(int z=0; z<a.rows(); z++)
                r(i,j) += a(i,z) * b(z,j);
    return r;
}

template <typename T>
Matrix<T> QRDescomposition<T>::scal_mat(Matrix<T> &a, int b){

    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            a(i,j) = a(i,j)*b;
    }
    return a;
}

#endif /* DESCOMPOSITION_QRDESCOMPOSITION_H_ */
