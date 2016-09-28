#ifndef _FFT_H_
#define _FFT_H_

#include <Eigen/Core>

using namespace Eigen;

void dft2d(const MatrixXd &m, MatrixXcd &f);
void dft2d(const MatrixXcd &m, MatrixXcd &f); 

void idft2d(const MatrixXcd &m, MatrixXd &f);
void idft2d(const MatrixXcd &m, MatrixXcd &f);

#endif
