#ifndef _DEBLUR_
#define _DEBLUR_

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;
using namespace std;

/*--------------------------------------------------------------------
 * class definition
 --------------------------------------------------------------------*/
class Deblur {

 public:
  Deblur(const char *filepath);
  virtual ~Deblur() {};

  ////  bool deblurImage();

  //// test
  void _test();

  Mat getImage(void) { return Image;};
  Mat getDebluredImage(void) { return debluredImage;};
  Mat getFilteredImage(void) { return filteredImage;};
  Mat getKernelImage(void) { return kernelImage;};


  ////  Mat getKernelImage(void) { return kernelImage;};

  void SharpFilter (const Mat &src, Mat &dst, const float k);
  void SharpFilter (const float k);

  void initKernelImage(const int size);

  void setParam(const float sigma_x, const float sigma_y, const float snr) {/*dummy*/};
  bool saveImage(const char *filepath) {return true; /*dummy*/};

  void setDebluredImage(Mat result);
  void setKernelImage(Mat kernel);
 private:
  Mat Image;
  Mat debluredImage;  
  Mat filteredImage;  
  Mat kernelImage;
  
  int height, width;
};
#endif
