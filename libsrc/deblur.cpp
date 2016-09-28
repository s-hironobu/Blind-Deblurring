/*
 * Deblur
 */

#include <iostream>
#include "deblur.h"

using namespace cv;
using namespace std;

enum {
  PANIC,
  ERR,
  WARNING
};


Deblur::Deblur(const char *filepath)
{
  try {
    Image = imread (filepath, CV_LOAD_IMAGE_COLOR);
    if (Image.empty ())
      throw ERR;

    height = Image.rows;
    width = Image.cols;
    debluredImage = Mat(height, width, CV_8U);
    filteredImage = Mat(height, width, CV_8U);

    debluredImage = Image.clone();
    filteredImage = Image.clone();    

    kernelImage = Mat::zeros(13,13,CV_8U);
    kernelImage.at<unsigned char>(6,6) = 255;

  }
  catch (int msg) {
    ;
  }
}

void Deblur::initKernelImage(const int size) {
  int center = size/2;
  kernelImage.release();
  kernelImage = Mat::zeros(size,size,CV_8U);
  kernelImage.at<unsigned char>(center,center) = 255;  
}

void Deblur::setDebluredImage(Mat result) {
    // TODO memory management
    if (result.depth() != 0 /* CV_8U */) {
      Mat tmp;
      result.convertTo(tmp, CV_8U, 255);
      tmp.copyTo(debluredImage);
    }
    else {
      result.copyTo(debluredImage);
    }
  };

void Deblur::setKernelImage(Mat kernel) {
    // TODO memory management
    if (kernel.depth() != 0 /* CV_8U */) {
      Mat tmp;
      kernel.convertTo(tmp, CV_8U, 255);
      tmp.copyTo(kernelImage);
    }
    else {
      kernel.copyTo(kernelImage);
    }
};

void Deblur::SharpFilter (const float k)
{
  SharpFilter(debluredImage, filteredImage, k);
}


void Deblur::SharpFilter (const Mat &src, Mat &dst, const float k)
{
  
  float KernelData[] = {
    -k/9.0f, -k/9.0f,           -k/9.0f,
    -k/9.0f, 1 + (8 * k)/9.0f,  -k/9.0f,
    -k/9.0f, -k/9.0f,           -k/9.0f,
  };
  Mat kernel = Mat (3, 3, CV_32F, KernelData);
  
  filter2D (src, dst, -1, kernel);
}
