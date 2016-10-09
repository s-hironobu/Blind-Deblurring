#ifndef _LIB_DEBLUR_
#define _LIB_DEBLUR_


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

using namespace cv;
using namespace std;
using namespace Eigen;


/*-------------------------------------------------------
 * enum
 -------------------------------------------------------*/
enum
  {
    DEBLUR_LOW,
    DEBLUR_MEDIUM,
    DEBLUR_HIGH
  };

/*-------------------------------------------------------
 * typedef
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
 * declare function
 -------------------------------------------------------*/
void dft2d(const MatrixXd &m, MatrixXcd &f);
void dft2d(const MatrixXcd &m, MatrixXcd &f); 
void idft2d(const MatrixXcd &m, MatrixXd &f);
void idft2d(const MatrixXcd &m, MatrixXcd &f);


/*--------------------------------------------------------------------
 * class definition
 --------------------------------------------------------------------*/
class Deblur {
 public:
  Deblur(const char *filepath);
  virtual ~Deblur() {};
  void _test();
  Mat getImage(void) { return Image;};
  Mat getDebluredImage(void) { return debluredImage;};
  Mat getFilteredImage(void) { return filteredImage;};
  Mat getKernelImage(void) { return kernelImage;};
  void SharpFilter (const Mat &src, Mat &dst, const float k);
  void SharpFilter (const float k);
  void initKernelImage(const int size);
  void setParam(const int mode);
  bool saveImage(const char *filepath) {return true; /*dummy*/};
  void setDebluredImage(Mat result);
  void setKernelImage(Mat kernel);
 private:
  Mat Image;
  Mat debluredImage;  
  Mat filteredImage;  
  Mat kernelImage;
    int height, width;
  int mode;
};

class BlindDeblur {
  public:
  BlindDeblur();
  ~BlindDeblur() {};

  virtual void setDeblur(Deblur *_deblur) { deblur = _deblur;};
  void exe(void);
  void exe(const Mat);
  void exe(const char *);

  void setMode(const int);
  void setPsfSize(const int);

  int getCount() {return count;};
  virtual bool checkInterrupt() { return false;};

  Mat debluredImage;
  Mat kernelImage;

  int count;

 protected:
  Deblur *deblur; 

private:
  param_t param;
  int mode;

  int sfix(const double);
  double sumabs2 (const MatrixXcd &);
  void imresize (const MatrixXd &, MatrixXd &, const int, const int);
  void imresize (const MatrixXd &, MatrixXd &, const float);
  void imresize (const MatrixXd &, MatrixXd &);
  void normalizeImage (MatrixXd &, double &, double &);
  void aLn (MatrixXd &, const MatrixXd &, const MatrixXd &, const vu_star_t);
  void asetLnorm (const double, double, double, vu_star_t *);
  void get_margins_mat (const MatrixXd &, int &, int &, int &, int &, const double);
  void make_mask_mat (const MatrixXd &, MatrixXd &, const double);
  void mat2gray (MatrixXd &);
  void bwmorph_clean (MatrixXd &);
  void set_Vh (MatrixXd &, const int, const int);
  void centerPSF (MatrixXd &, const param_t);

  void PSFestimaLnoRgrad (MatrixXd &, MatrixXcd &, param_t, const int);
  void Ustep (const double, MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXcd &,
	      MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &,
	      MatrixXd &, MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &,
	      const int, const int, param_t, double, const int, const int);
  void Hstep (double, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXd &, MatrixXd &, MatrixXd &,
	      MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXd &, const int, const int, const int, const int, param_t, double, const int, const int);
  void set_vrange (const MatrixXd &, vrange_t &);
  void uConstr (MatrixXcd &, const vrange_t);
  void fftCGSRaL (MatrixXd &, MatrixXd &, MatrixXcd &, const param_t);
  MatrixXd getROI (const MatrixXd &, const param_t);
  void createROI(MatrixXd &, const int, const param_t);
  void doublePSF(MatrixXd &);
  void divide_rgb_images(Mat, MatrixXd &, MatrixXd &, MatrixXd &,
			 double &, double &, double &, double &, double &, double &);
  void make_image(Mat &, const MatrixXcd, const MatrixXcd, const MatrixXcd, const double, 
		  const double, const double, const double, const double, const double);
  void show_psf(const MatrixXd &);
  void matNormalize(MatrixXd &);
  void sharp_filter(Mat &, const param_t);
  void init_param (param_t *);
  void init_param ();  

  // workload
  int totalWorkload;
  int weight;
  int counter;
  void calcTotalWorkload(param_t *);

  void countUp(void);
  void changeWeight(const int);
};

#endif
