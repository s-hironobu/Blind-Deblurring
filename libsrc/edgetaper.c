/*-------------------------------------------------------
 * edgetaper
 *
 *
 -------------------------------------------------------*/
#include "opencv2/core/core.hpp"
#include <Eigen/Core>

using namespace Eigen;

static void
create_gaussian_kernel(MatrixXd & psf, double sigma)
{
  int psf_x = psf.rows();
  int psf_y = psf.cols();

  double factor = 1.0 / (sqrt(2.0 * M_PI) * sigma);
  double hh = psf_x / 2.0;
  double hw = psf_y / 2.0;
  
  double tmp_x, tmp_y;
  double sum = 0.0;
  
  for (int x = 0; x < psf_x; x++) {
    tmp_x = (x - hh) / sigma;
    for (int y = 0; y < psf_y; y++) {
      tmp_y = (y - hw) / sigma;
      psf(x,y) = exp(-(tmp_x * tmp_x + tmp_y * tmp_y) / 2.0) * factor;
      sum += psf(x,y);
    }
  }
  if (sum > 0.0)
    psf.array()  /= sum;
}

static void
normalized_autocorrelation(const MatrixXd psf, MatrixXd & xpsf) {

  int psf_x = (int)psf.rows();
  int psf_y = (int)psf.cols();
  int xpsf_x = (int)xpsf.rows();
  int xpsf_y = (int)xpsf.cols();

  assert((2*psf_x -1) == xpsf_x);
  assert((2*psf_y -1) == xpsf_x);

  //      ipsf(x-u,y-v) = psf(psf_x -1 - x + u, psf_y - 1 - y + v);
  /*
   * C(x,y) = Sigma_u Sigma_v A(u,v) * B(x-u+1,y-v+1)    
   */
  double max = DBL_MIN;
  for (int x = 0; x < xpsf_x; x++)
    for (int y = 0; y < xpsf_y; y++) { 
      double c = 0.0;
      for (int u = 0; u < psf_x; u++) {
	for (int v = 0; v < psf_y; v++) {
	  if (0<=(x-u) && 0<=(y-v)
	      && (x-u)<psf_x && (y-v)<psf_y) {
	    //c += psf(u,v) * ipsf(x-u,y-v);
	    c += psf(u,v) * psf(psf_x -1 - x + u, psf_y - 1 - y + v);
	  }
	}
      }
      xpsf(x,y) = c;
      if (max < c)	max = c;      
    }
  /*
   * normarize
   */
  if (max != 0)
    xpsf.array() /= max;
}

static void
create_weights(MatrixXd &weights, const MatrixXd psf)
{
  int psf_x = (int)psf.rows();
  int psf_y = (int)psf.cols();
  int w_x = (int)weights.rows();
  int w_y = (int)weights.cols();

  assert(2*psf_x < w_x);
  assert(2*psf_y < w_y);

  int mergin_x = psf_x - 1;
  int mergin_y = psf_y - 1;

  MatrixXd xpsf(2*psf_y -1, 2*psf_y -1);
  normalized_autocorrelation(psf, xpsf);

  for (int x = 0; x < w_x; x++)
    for (int y = 0; y < w_y; y++) {
      //
      if (x < mergin_x)
	if (y < mergin_y)
	  weights(x,y) = xpsf(x,y); // (I)
	else if (mergin_y <= y && y < (w_y - mergin_y))
	  weights(x,y) = xpsf(x, mergin_y); // (IV)
	else
	  weights(x,y) = xpsf(x, y - (w_y - 2*mergin_y)+1); // (VII)
      //
      else if (mergin_x <= x && x < (w_x - mergin_x))
	if (y < mergin_y)
	  weights(x,y) = xpsf(mergin_x, y); // (II)
	else if (mergin_y <= y && y < (w_y - mergin_y))
	  weights(x,y) = 1.0; //  (V)
	else
	  weights(x,y) = xpsf(mergin_x , y - (w_y - 2*mergin_y)+1); // (IIX)
	//
      else
	if (y < mergin_y)
	  weights(x,y) = xpsf(x - (w_x - 2*mergin_x) +1, y); // (III)
	else if (mergin_y <= y && y < (w_y - mergin_y))
	  weights(x,y) = xpsf(x - (w_x - 2*mergin_x) +1, mergin_y); // (VI) 
	else
	  weights(x,y) = xpsf(x - (w_x - 2*mergin_x)+1, y - (w_y - 2*mergin_y)+1); // (IX)
    }  
}

static void
create_blurred (MatrixXcd & blurred, const MatrixXcd img, const MatrixXd psf)
{
  int img_x = (int)img.rows();
  int img_y = (int)img.cols();
  int psf_x = (int)psf.rows();
  int psf_y = (int)psf.cols();
  
  int m_x = (((int)psf_x+1) >> 1) -1;
  int m_y = (((int)psf_y+1) >> 1) -1;

  for (int x = 0; x < img_x; x++)
    for (int y = 0; y < img_y; y++) { 

      if (((psf_x < x && x <= (img_x - psf_x))  && ((psf_y < y && y <= (img_y - psf_y))))) {
	blurred(x,y) = img(x,y);
      }
      else {
	for (int i = 0; i < psf_x; i++) {
	  for (int j = 0; j < psf_y; j++) {
	    int u = i - m_x;
	    int v = j - m_y;
	    if ((0<=(x+u)) && (0<=(y+v))
		&& ((x+u) < img_x) && ((y+v) < img_y)) {
	      blurred(x,y) += psf(i,j) * img(x+u,y+v);
	    }
	  }
	}      
      }      
    }
}

void
edgetaper(MatrixXcd & edge, const MatrixXcd img)
{
  int edge_x = (int)edge.rows();
  int edge_y = (int)edge.cols();
  int img_x = (int)img.rows();
  int img_y = (int)img.cols();

  int psf_x, psf_y;
  psf_x = psf_y = 10;

  MatrixXd psf(psf_x, psf_y);
  create_gaussian_kernel(psf, (float)1.0);

  assert(edge_x == img_x && edge_y == img_y);
  assert(2*psf_x < edge_x && 2*psf_y < edge_y);

  /* blurred
   */
  MatrixXcd blurred(img_x, img_y);
  create_blurred(blurred, img, psf);

  /* weights
   */
  MatrixXd weights(img_x, img_y);
  create_weights(weights, psf);

  /*
   */
  for (int x = 0; x < img_x; x++) {
    for (int y = 0; y < img_y; y++) {
      if (weights(x,y) != 1.0)
	edge(x,y) = img(x,y)*weights(x,y) + blurred(x,y) *((double)1.0 - weights(x,y));
      else
	edge(x,y) = img(x,y);
    }
  }
}

