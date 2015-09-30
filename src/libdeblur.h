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
    DEBLUR_HIGH,
    DEBLUR_MEDIUM,
    DEBLUR_LOW
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
  int base_8c9afa73d5fb2e5a964ee212c8ae11eb;
  int base_1b7de9732286185efd1dbce52a1f29ba;
  int base_9d2e966f4abf99089fa91c874cbbc4aa;
  int base_c8bb87c35d5c8016039eb2c63f5055a1;
  int base_c8672b53175c78bd40dfa26a7aae49ec;
  double base_05b048d7242cb7b8b57cfa3b1d65ecea;
  double base_83aa590baf7250511fba793140bd1ea0;
  double base_2ec348ba37e0003b20442504ceb32648;
  double base_3643e51d870f4849a9dfe990c779ad9d;
  double base_48614ae9d92989d76e508df68f7399b6;
  double base_8c4128e380caf5f0e88a870e04e8c397;
  double base_5f17275558315808e826bd4d916f9917;
  double base_90c0215603f7090fd20ba9e4cfecb70f;
  double base_085c6920b37fd444ee73a0ebe8e929c6;
  int base_331cba6400bc264dc96dd7549b144eb7;
  int base_91155c5fcc54f25919f9ac1461ca7299;
  int base_c1c7c2b22243f2eaea86834e3d7bb451;
  int base_4508a9419f4e73e8274f03db1a5b333f;
  double base_8baf4edfac39cf03451a836f1fd0a079;
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
  void base_1872614a65a70d51910ee7ced46d9bf5 (const MatrixXd &, MatrixXd &, const int, const int);
  void base_1872614a65a70d51910ee7ced46d9bf5 (const MatrixXd &, MatrixXd &, const float);
  void base_1872614a65a70d51910ee7ced46d9bf5 (const MatrixXd &, MatrixXd &);
  void base_5b18cd4b4b96f7299e4097fda90288e3 (MatrixXd &, double &, double &);
  void base_2fea1f166d05e7dde8763061bb16e801 (MatrixXd &, const MatrixXd &, const MatrixXd &, const vu_star_t);
  void base_85745d479dd5e1cd108b054d8a48620b (const double, double, double, vu_star_t *);
  void base_7cccd49cb2e59a1d0632bd4befb14c67 (const MatrixXd &, int &, int &, int &, int &, const double);
  void base_edfcc099f60103785bd7668cdbb49cca (const MatrixXd &, MatrixXd &, const double);
  void base_8a42df5b19795bad5870e5a2cacec755 (MatrixXd &);
  void base_cbb2ace6efa94a7fe34ca9ee1b542e53 (MatrixXd &);
  void base_6667ce9bd69d8446067e41f669194952 (MatrixXd &, const int, const int);
  void base_6de88db41628e1fd44e929ffcb14070f (MatrixXd &, const param_t);

  void PSFestimbase_2fea1f166d05e7dde8763061bb16e801oRgrad (MatrixXd &, MatrixXcd &, param_t, const int);
  void base_2079a4b8e2329f166ff90bc9604b1956 (const double, MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXcd &,
	      MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &,
	      MatrixXd &, MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXd &, MatrixXd &,
	      const int, const int, param_t, double, const int, const int);
  void base_4e59fbcfff3e3dd649cdf148726bffcb (double, MatrixXd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXd &, MatrixXd &, MatrixXcd &, MatrixXd &, MatrixXd &, MatrixXd &,
	      MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &, MatrixXcd &,
	      MatrixXd &, const int, const int, const int, const int, param_t, double, const int, const int);
  void base_3c005cc3ace054db49702f3003191eee (const MatrixXd &, vrange_t &);
  void base_00adb5c87e5ac2ab7ac7c442d1ccc95a (MatrixXcd &, const vrange_t);
  void base_930ca8b114a91ebda3e9e9697d00cb48 (MatrixXd &, MatrixXd &, MatrixXcd &, const param_t);
  MatrixXd base_c6dbff44d02a742559d452259555580c (const MatrixXd &, const param_t);
  void base_a813d46a372ef1707696d5afa5938ac8(MatrixXd &, const int, const param_t);
  void base_dbd3d60ff27d1ece1c15f497380749f2(MatrixXd &);
  void base_d16862fd866aec72be2ef18eba71e2a0(Mat, MatrixXd &, MatrixXd &, MatrixXd &,
			 double &, double &, double &, double &, double &, double &);
  void base_a615dbaa3cce1eb9f55ff8fd4c312019(Mat &, const MatrixXcd, const MatrixXcd, const MatrixXcd, const double, 
		  const double, const double, const double, const double, const double);
  void base_2dd381f8cfd647e52417b13c61c06e7a(const MatrixXd &);
  void base_1b77cbdb44bc9c713e38d4a055e1247a(MatrixXd &);
  void base_bc63da0ff03c844bf2b48f01a9e14aa5(Mat &, const param_t);
  void base_841040e67396f60862af8530508993a8 (param_t *);
  void base_841040e67396f60862af8530508993a8 ();  

  // workload
  int base_cc7bcf5e6a4856507f072bee6f077389;
  int base_7edabf994b76a00cbc60c95af337db8f;
  int base_886bb73b3156b0aa24aac99d2de0b238;
  void base_ce2698d4f64659528f0613662206f59d(param_t *);

  void base_dbd14b23cdaf2dc59b90da4493314837(void);
  void base_7ba3b0eea66eec11b838f65475363d1b(const int);
};

#endif
