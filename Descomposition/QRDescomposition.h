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
    bool solveQR(const Matrix<T> &A, vector<T> &x, const vector<T> &b);
private:
    Matrix<T> ext_prod(Matrix<T> &a, Matrix<T> &b);
    Matrix<T> mat_prod(Matrix<T> &a, Matrix<T> &b);
    Matrix<T> scal_mat(Matrix<T> &a, int b);
    Matrix<T> transpose(Matrix<T> &a);
    void printMatrix(anpi::Matrix<T> &m);
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
        printMatrix(P);
        cout<<endl;
        printMatrix(R);
        cout<<endl;
        Matrix<T> R1=R;
        R = mat_prod(P,R);
        printMatrix(R);
        R1=P*R1;
         cout<<endl;
        printMatrix(R1);
         cout<<endl;
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

template <typename T>
Matrix<T> QRDescomposition<T>::transpose(Matrix<T> &a){
    T x = 0;
    Matrix<T> r(a.cols(),a.rows(),x);
    for(int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            r(j,i) = a(i,j);
    }
    return r;
}


template <typename T>
bool QRDescomposition<T>::solveQR(const Matrix<T> &A, vector<T> &x, const vector<T> &b){
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

template <typename T>
void QRDescomposition<T>::printMatrix(anpi::Matrix<T> &m){
    for(int i = 0; i < m.rows(); i++){
        cout << "|\t";
        for(int j = 0; j < m.cols(); j++)
            cout << "[" << m[i][j] << "]\t";
        cout << "|" <<endl;
    }
}

#endif /* DESCOMPOSITION_QRDESCOMPOSITION_H_ */
