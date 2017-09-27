/*
 * MatrixException.hpp
 *
 *  Created on: 25 de set. de 2017
 *      Author: kevin
 */

#ifndef MATRIX_MATRIXEXCEPTION_HPP_
#define MATRIX_MATRIXEXCEPTION_HPP_



#include <iostream>
#include <stdexcept>

class MatrixException : public std::exception
{
public:

	explicit MatrixException(){}

	virtual const char* what() const throw(){
		return "Number of columns of the first matrix is different to the number of rows of the second";
	}
};


#endif /* MATRIX_MATRIXEXCEPTION_HPP_ */
