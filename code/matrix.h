//Matrix class
#ifndef __MATRIX_H__
#define __MATRIX_H__

//#include <math>

class matrix {
	public:
	double a0,a1,a2,a3;
	matrix(double b0, double b1, double b2, double b3);
	matrix(const matrix& M);
	matrix();
	// Trace
	double tr();
	// inverse matrix
	void inverse();
	matrix conj();
	//Projection to SU(2)
	void proj();
};

//overloaded matrix operations
matrix operator+ (matrix A, matrix B);
matrix operator- (matrix A, matrix B);
matrix operator* (double x, matrix A);
matrix operator* (matrix A, double x);
matrix operator* (matrix A, matrix B);
#endif
