#ifndef _BLIND_DEBLUR_
#define _BLIND_DEBLUR_

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdio.h>
#include <math.h>
#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "fft.h"

#include "deblur.h"

using namespace cv;
using namespace std;
using namespace Eigen;

/*-------------------------------------------------------
 * typedef
 *
 *
 -------------------------------------------------------*/

typedef struct vrange
{
  double min;
  double max;
} vrange_t;

typedef struct param
{
  int mode;

  int psf_size_x;
  int psf_size_y;

  int maxROIsize_r;
  int maxROIsize_c;

  int MSlevels;

  double gamma;
  double Lp;

  double beta_h;
  double alpha_h;
  double centering_threshold;

  double beta_u;
  double alpha_u;

  double gamma_nonblind;
  double beta_u_nonblind;
  int Lp_nonblind;

  int maxiter_u;
  int maxiter_h;
  int maxiter;
  double ccreltol;

  double filter;

} param_t;

typedef struct vu_star
{
  double v_star;
  double u_star;
} vu_star_t;


typedef struct setLnorm
{
  double x;
  double alpha;
  double beta;
  double q;
} setLnorm_t;


/*-------------------------------------------------------
 * enum
 *
 *
 -------------------------------------------------------*/
enum
  {
    DEBLUR_HIGH,
    DEBLUR_MEDIUM,
    DEBLUR_LOW
  };


/*-------------------------------------------------------
 * class
 *
 *
 -------------------------------------------------------*/
class BlindDeblur {
  
public:
  BlindDeblur();
  ~BlindDeblur() {};

  virtual void setDeblur(Deblur *_deblur) { deblur = _deblur;};
  void exe(void);
  void exe(const Mat img);
  void exe(const char *filename);

  void setMode(const int mode);
  void setPsfSize(const int size);

  int getCount() {return count;};
  virtual bool checkInterrupt() { return false;};

  Mat debluredImage;
  Mat kernelImage;

  int count;

 protected:
  Deblur *deblur; 

private:
  param_t param;

  int sfix(const double d);
  double sumabs2 (const MatrixXcd & m);
  void imresize (const MatrixXd &in, MatrixXd &out, const int sx, const int sy);
  void imresize (const MatrixXd &in, MatrixXd &out, const float m);
  void imresize (const MatrixXd &in, MatrixXd &out);
  void normalizeImage (MatrixXd & m, double &min, double &max);
  void aLn (MatrixXd & V, const MatrixXd & DU, const MatrixXd & normDU, const vu_star_t r);
  void asetLnorm (const double q, double alpha, double beta, vu_star_t * r);
  void get_margins_mat (const MatrixXd & m, int &m_left, int &m_right, int &m_top,
			int &m_bottom, const double threshold);
  void make_mask_mat (const MatrixXd & m, MatrixXd & mask, const double threshold);
  void mat2gray (MatrixXd & m);
  void bwmorph_clean (MatrixXd & mask);
  void set_Vh (MatrixXd & vh, const int rows, const int cols);
  void centerPSF (MatrixXd & H, const param_t param);
  void PSFestimaLnoRgrad (MatrixXd & h, MatrixXcd & ROI, param_t param, const int L);
  void Ustep (const double q,
	      MatrixXd & H, MatrixXd & U, MatrixXcd & FeGu,
	      MatrixXcd & FU, MatrixXcd & FUx, MatrixXcd & FUy,
	      MatrixXcd & FDx, MatrixXcd & FDy,
	      MatrixXd & Vx, MatrixXd & Vy,
	      MatrixXd & Bx, MatrixXd & By,
	      MatrixXcd & DTD,
	      MatrixXcd & FHS, MatrixXcd & FHTH, MatrixXcd & FGs,
	      MatrixXcd & FUp, MatrixXcd & b,
	      MatrixXd & xD, MatrixXd & yD,
	      MatrixXd & xDm, MatrixXd & yDm, MatrixXd & nDm,
	      MatrixXcd & tmp1, MatrixXcd & tmp2, 
	      MatrixXd & tmpr1, MatrixXd & tmpr2,
	      const int usize_r, const int usize_c, param_t param, double gamma, const int L, const int ml);
  void Hstep (double q,
	      MatrixXd & H,
	      MatrixXcd & FeGx, MatrixXcd & FeGy,
	      MatrixXcd & FUx, MatrixXcd & FUy,
	      MatrixXd & Vh, MatrixXd & Bh,
	      MatrixXcd & FHp, MatrixXd & hI,
	      MatrixXd & hIm, MatrixXd & nIm,
	      MatrixXcd & FUD, MatrixXcd & FUTU, MatrixXcd & FH,
	      MatrixXcd & b, MatrixXcd & tmp, MatrixXcd & tmp1,
	      MatrixXd & tmpr1,
	      const int usize_r, const int usize_c,
	      const int hsize_r, const int hsize_c, 
	      param_t param, double gamma, const int L, const int ml);
  void set_vrange (const MatrixXd & m, vrange_t & vrange);
  void uConstr (MatrixXcd & m, const vrange_t vrange);
  void fftCGSRaL (MatrixXd & G, MatrixXd & H, MatrixXcd & U, const param_t param);
  MatrixXd getROI (const MatrixXd & m, const param_t param);
  void createROI(MatrixXd &tmp, const int L, const param_t param);
  void doublePSF(MatrixXd & h);
  void divide_rgb_images(Mat color_img,
			 MatrixXd & r_src, MatrixXd & g_src, MatrixXd & b_src,
			 double & r_min, double & r_max,
			 double & g_min, double & g_max,
			 double & b_min, double & b_max);
  void make_image(Mat & result, const MatrixXcd rU, const MatrixXcd gU, const MatrixXcd bU, 
		  const double r_min, const double r_max,
		  const double g_min, const double g_max,
		  const double b_min, const double b_max);
  void show_psf(const MatrixXd &h);
  void matNormalize(MatrixXd & h);
  void sharp_filter(Mat &result, const param_t param);
  void init_param (param_t * param);
  void init_param ();  

  // workload
  int totalWorkload;
  int weight;
  int counter;
  void calcTotalWorkload(param_t * param);
  void countUp(void);
  void changeWeight(const int weight);
};

#endif
