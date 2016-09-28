#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include "fftw3.h"

using namespace std;
using namespace Eigen;

#define _FFTW3

#ifdef _FFTW3
void 
dft2d(const MatrixXd &m, MatrixXcd &f) 
{
  fftw_complex *in  = NULL;
  fftw_complex *out = NULL;

  fftw_plan p       = NULL;

  int idx;
  int size = m.rows() * m.cols();
  
  size_t mem_size = sizeof(fftw_complex) * size;
  in  = (fftw_complex*)fftw_malloc( mem_size );
  out = (fftw_complex*)fftw_malloc( mem_size );
  
  if( !in || !out ){
    fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
    exit(-1);
  }
  
  p = fftw_plan_dft_2d(m.rows(), m.cols(), in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  for(int j = 0; j < m.rows(); j++) {
    for(int i = 0; i < m.cols(); i++) {
      idx = m.cols() * j + i; // column-major alignment
      in[idx][0] = m(j, i);
      in[idx][1] = 0.0;
    }
  }
  
  fftw_execute(p);

  for (int j = 0; j < m.rows(); j++) {
    for (int i = 0; i < m.cols(); i++){
      idx = m.cols() * j + i;
      f(j,i) = complex<double>(out[idx][0], out[idx][1]);
    }
  }
  if( p   ) fftw_destroy_plan(p);
  if( in  ) fftw_free(in);
  if( out ) fftw_free(out);


}

void 
dft2d(const MatrixXcd &m, MatrixXcd &f) 
{
  fftw_complex *in  = NULL;
  fftw_complex *out = NULL;

  fftw_plan p       = NULL;

  int idx;
  int size = m.rows() * m.cols();
  
  size_t mem_size = sizeof(fftw_complex) * size;
  in  = (fftw_complex*)fftw_malloc( mem_size );
  out = (fftw_complex*)fftw_malloc( mem_size );
  
  if( !in || !out ){
    fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
    exit(-1);
  }
  
  p = fftw_plan_dft_2d(m.rows(), m.cols(), in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  for(int j = 0; j < m.rows(); j++) {
    for(int i = 0; i < m.cols(); i++) {
      idx = i + m.cols() * j; // column-major alignment
      in[idx][0] = real(m(j, i));
      in[idx][1] = imag(m(j, i));;
    }
  }
  
  fftw_execute(p);

  for (int j = 0; j < m.rows(); j++){
    for (int i = 0; i < m.cols(); i++) {
      idx = i +  m.cols() * j;
      f(j,i) = complex<double>(out[idx][0], out[idx][1]);
    }
  }
  if( p   ) fftw_destroy_plan(p);
  if( in  ) fftw_free(in);
  if( out ) fftw_free(out);
}

#else

void
dft2d (const MatrixXd &m, MatrixXcd &f)
{
  Mat tmp;
  Mat M = Mat(m.rows (), m.cols (), CV_64F);
  Mat planes[2] = { Mat_ < double >(M), Mat::zeros (M.size (), CV_64F) };
  
  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      {
	planes[0].at < double >(x, y) = (double)m (x, y);
	planes[1].at < double >(x, y) = (double)0.0;
      }
  
  merge (planes, 2, tmp);
  dft (tmp, tmp, DFT_COMPLEX_OUTPUT);
  
  vector < cv::Mat_<double> > F;
  Mat f_r = Mat(m.rows (), m.cols (), CV_64F);
  Mat f_i = Mat(m.rows (), m.cols (), CV_64F);

  split (tmp, F);
  f_r = F[0];
  f_i = F[1];
  
  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++) {
      f (x, y) =
	complex < double >(f_r.at < double >(x, y), f_i.at < double >(x, y));
    }
}

void
dft2d (const MatrixXcd & m, MatrixXcd & f)
{
  Mat tmp;
  Mat M = Mat(m.rows (), m.cols (), CV_64F);
  Mat planes[2] = { Mat_ < double >(M), Mat::zeros (M.size (), CV_64F) };
  
  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      {
	planes[0].at < double >(x, y) = (double)real (m (x, y));
	planes[1].at < double >(x, y) = (double)imag (m (x, y));
      }
  
  merge (planes, 2, tmp);
  dft (tmp, tmp, DFT_COMPLEX_OUTPUT);
  
  vector < cv::Mat_<double> > F;
  Mat f_r = Mat(m.rows (), m.cols (), CV_64F);
  Mat f_i = Mat(m.rows (), m.cols (), CV_64F);

  split (tmp, F);
  f_r = F[0];
  f_i = F[1];
  
  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++) {
      f (x, y) =
	complex < double >(f_r.at < double >(x, y), f_i.at < double >(x, y));
    }
}
#endif

#ifdef _FFTW3

void 
idft2d(const MatrixXcd &m, MatrixXd &f) 
{
  fftw_complex *in  = NULL;
  fftw_complex *out = NULL;

  fftw_plan p       = NULL;

  int idx;
  int size = m.rows() * m.cols();
  
  size_t mem_size = sizeof(fftw_complex) * size;
  in  = (fftw_complex*)fftw_malloc( mem_size );
  out = (fftw_complex*)fftw_malloc( mem_size );
  
  if( !in || !out ){
    fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
    exit(-1);
  }
  
  // !! row-major alignment is recommended, but here, column-major.
  p = fftw_plan_dft_2d(m.rows(), m.cols(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  for(int j = 0; j < m.rows(); j++) {
    for(int i = 0; i < m.cols(); i++) {
      idx = m.cols() * j + i; // column-major alignment
      in[idx][0] = real(m(j, i));
      in[idx][1] = imag(m(j, i));
    }
  }
  
  fftw_execute(p);
  
  // output is DC exchanged and scaled.
  double scale = 1. / (m.rows() * m.cols());
  for (int j = 0; j < m.rows(); j++) {
    for (int i = 0; i < m.cols(); i++){
      idx = m.cols() * j + i;
      f(j,i) = out[idx][0] * scale;
    }
  }
  
  if( p   ) fftw_destroy_plan(p);
  if( in  ) fftw_free(in);
  if( out ) fftw_free(out);
}

void 
idft2d(const MatrixXcd &m, MatrixXcd &f) 
{
  fftw_complex *in  = NULL;
  fftw_complex *out = NULL;

  fftw_plan p       = NULL;

  int idx;
  int size = m.rows() * m.cols();
  
  size_t mem_size = sizeof(fftw_complex) * size;
  in  = (fftw_complex*)fftw_malloc( mem_size );
  out = (fftw_complex*)fftw_malloc( mem_size );
  
  if( !in || !out ){
    fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
    exit(-1);
  }
  
  // !! row-major alignment is recommended, but here, column-major.
  p = fftw_plan_dft_2d(m.rows(), m.cols(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  for(int j = 0; j < m.rows(); j++) {
    for(int i = 0; i < m.cols(); i++) {
      idx = m.cols() * j + i; // column-major alignment
      in[idx][0] = real(m(j, i));
      in[idx][1] = imag(m(j, i));
    }
  }
  
  fftw_execute(p);
  
  // output is DC exchanged and scaled.
  double scale = 1. / (m.rows() * m.cols());
  for (int j = 0; j < m.rows(); j++) {
    for (int i = 0; i < m.cols(); i++){
      idx = m.cols() * j + i;
      f(j,i) = complex<double>(out[idx][0] * scale, out[idx][1] * scale);
    }
  }
  
  if( p   ) fftw_destroy_plan(p);
  if( in  ) fftw_free(in);
  if( out ) fftw_free(out);
}

#else

void
idft2d (const MatrixXcd & m, MatrixXd & f)
{
  Mat tmp;
  Mat M (m.rows (), m.cols (), CV_64F);
  Mat planes[] = { Mat_ < double >(M), Mat::zeros (M.size (), CV_64F) };

  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      {
	planes[0].at < double >(x, y) = real (m (x, y));
	planes[1].at < double >(x, y) = imag (m (x, y));
      }
  merge (planes, 2, tmp);
  idft (tmp, tmp, DFT_SCALE);

  vector < cv::Mat > F;
  Mat f_r, f_i;
  split (tmp, F);
  f_r = F[0];
  f_i = F[1];

  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      f (x, y) = f_r.at < double >(x, y);
}

void
idft2d (const MatrixXcd & m, MatrixXcd & f)
{

  Mat tmp;
  Mat M (m.rows (), m.cols (), CV_64F);
  Mat planes[] = { Mat_ < double >(M), Mat::zeros (M.size (), CV_64F) };

  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      {
	planes[0].at < double >(x, y) = real (m (x, y));
	planes[1].at < double >(x, y) = imag (m (x, y));
      }
  merge (planes, 2, tmp);
  idft (tmp, tmp, DFT_SCALE);

  vector < cv::Mat > F;
  Mat f_r, f_i;
  split (tmp, F);
  f_r = F[0];
  f_i = F[1];

  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      f (x, y) =
	complex < double >(f_r.at < double >(x, y), f_i.at < double >(x, y));
}

#endif
