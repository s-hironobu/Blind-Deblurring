#include <Eigen/Core>

using namespace Eigen;

MatrixXd create_mat (const int rows, const int cols)
{
  return MatrixXd::Zero (rows, cols);
}

MatrixXcd create_complex_mat (const int rows, const int cols)
{
  return MatrixXcd::Zero (rows, cols);
}

void copy_cmat_2_mat (const MatrixXcd & m, MatrixXd & h)
{
  assert(m.rows() == h.rows() && m.cols() == h.cols());
  h = m.block(0, 0, m.rows(), m.cols()).real();
}

void copy_mat_2_cmat (const MatrixXd & a, MatrixXcd & b, const int rows, const int cols)
{
  // fftCGSRaL (G)
  b.setZero();
  assert(rows <= b.rows() && cols <= b.cols());
  b.block(0,0,rows-1,cols-1).real().array() = a.block(0,0,rows-1,cols-1).array();
}

void copy_mat_2_mat_zeros (const MatrixXd & a, MatrixXd & b, const int rows, const int cols)
{
  // Ustep (H)
  b.setZero();
  assert(rows <= b.rows() && cols <= b.cols());
  b.block(0,0,rows-1,cols-1).array() = a.block(0,0,rows-1,cols-1).array();
}


void copy_mat_2_cmat_zeros (const MatrixXd & a, MatrixXcd & b, const int rows, const int cols)
{
  // Hstep (H), fftCGSRaL(H)
  b.setZero();
  assert(rows <= b.rows() && cols <= b.cols());
  b.block(0,0,rows-1,cols-1).real().array() = a.block(0,0,rows-1,cols-1).array();
}


void copy_cmat_2_cmat_zeros (const MatrixXcd & a, MatrixXcd & b, const int rows, const int cols)
{
  // fftCGSRaL (H)
  b.setZero();
  assert(rows <= b.rows() && cols <= b.cols());
  b.block(0,0,rows-1,cols-1).array() = a.block(0,0,rows-1,cols-1).array();
}
