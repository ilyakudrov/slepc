//Matrix class

#include <cmath>
#include "matrix.h"

matrix::matrix(double b0, double b1, double b2, double b3) { a0=b0; a1=b1; a2=b2; a3=b3; }
matrix::matrix(const matrix& M) { a0=M.a0; a1=M.a1; a2=M.a2; a3=M.a3; }
matrix::matrix() {a0 = 1; a1 = 0; a2 = 0; a3 = 0;}

double matrix::tr() {return 2*a0; }
void matrix::inverse() { 
	double rho=a0*a0+a1*a1+a2*a2+a3*a3;
	a0=a0/rho; a1=-a1/rho; a2=-a2/rho; a3=-a3/rho; 
}
matrix matrix::conj() {
	return matrix(a0, -a1, -a2, -a3);
}
void matrix::proj(){
	double rho = a0*a0+a1*a1+a2*a2+a3*a3;
	a0 = a0/powf(rho, 0.5); a1 = a1/powf(rho, 0.5); a2 = a2/powf(rho, 0.5); a3 = a3/powf(rho, 0.5);
}

matrix operator+ (matrix A, matrix B) {
	return matrix(A.a0+B.a0,A.a1+B.a1,A.a2+B.a2,A.a3+B.a3); };
matrix operator- (matrix A, matrix B) {
	return matrix(A.a0-B.a0,A.a1-B.a1,A.a2-B.a2,A.a3-B.a3); };
matrix operator* (double x, matrix A) {
	return matrix(A.a0*x,A.a1*x,A.a2*x,A.a3*x); };
matrix operator* (matrix A, double x) {
	return matrix(A.a0*x,A.a1*x,A.a2*x,A.a3*x); };
matrix operator* (matrix A, matrix B) {
	return matrix(A.a0*B.a0-A.a1*B.a1-A.a2*B.a2-A.a3*B.a3,
		A.a0*B.a1+B.a0*A.a1+A.a3*B.a2-A.a2*B.a3,
		A.a0*B.a2+B.a0*A.a2+A.a1*B.a3-A.a3*B.a1,
		A.a0*B.a3+B.a0*A.a3+A.a2*B.a1-A.a1*B.a2); };
