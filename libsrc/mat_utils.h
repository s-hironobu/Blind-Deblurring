#ifndef _MAT_UTILS_H_
#define _MAT_UTILS_H_

#include <Eigen/Core>

using namespace Eigen;

MatrixXd create_mat (const int rows, const int cols);
MatrixXcd create_complex_mat (const int rows, const int cols);

void copy_cmat_2_mat (const MatrixXcd & m, MatrixXd & h);
void copy_mat_2_cmat (const MatrixXd & a, MatrixXcd & b, const int rows, const int cols);
void copy_mat_2_mat_zeros (const MatrixXd & a, MatrixXd & b, const int rows, const int cols);
void copy_mat_2_cmat_zeros (const MatrixXd & a, MatrixXcd & b, const int rows, const int cols);
void copy_cmat_2_cmat_zeros (const MatrixXcd & a, MatrixXcd & b, const int rows, const int cols);

#endif
