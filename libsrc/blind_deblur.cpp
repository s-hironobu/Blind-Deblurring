/*
 * blind_deblur <image>
 */
#include "blind_deblur.h"

#include "mat_utils.h"
#include "edgetaper.h"

#define _TIME_

#ifdef _TIME_
#include <sys/time.h>

static struct timeval stat_data_begin;
static struct timeval stat_data;

static double get_interval(struct timeval bt, struct timeval et)
{
  double b, e;
  b = bt.tv_sec + (double) bt.tv_usec * 1e-6;
  e = et.tv_sec + (double) et.tv_usec * 1e-6;
  return e - b;
}

static void START_TIME()
{
  gettimeofday(&stat_data_begin, NULL);
}

static void SHOW_TIME(const char *s, const int i)
{
  gettimeofday(&stat_data, NULL);
  printf("DIFF=%e(sec)  %s(%d)\n", get_interval(stat_data_begin, stat_data), s, i);
  fflush(stdout);
}

#endif // _TIME_


#define smax(a,b) \
  ((a>b)? a:b)
#define smin(a,b) \
  ((a>b)? b:a)

#define _EDGETAPER

//#define _VERBOSE0

#define MAX_LEVEL 4

/*-------------------------------------------------------
 * enum
 *
 *
 -------------------------------------------------------*/
enum {
  ERR_1,  // file not found
  ERR_2,
  ERR_3,
  ERR_4,

  INTR_1,
  INTR_2,
  INTR_3,
  INTR_4
};

/*-------------------------------------------------------
 * constructor
 *
 *
 *-------------------------------------------------------*/
BlindDeblur::BlindDeblur() { 
  init_param();
  deblur = NULL;
};


/*-------------------------------------------------------
 * public method
 *
 *
 *-------------------------------------------------------*/
void BlindDeblur::exe(const char *filename)
{
  //  param_t param;
  //  init_param (&param);

  try {
    Mat img = imread(filename, CV_LOAD_IMAGE_COLOR);
    if (img.empty ())
      throw ERR_1;
    
    exe(img);
    
  }catch(int err) {
    ;;//    throw err;
  }

#ifdef _TEST_  
#ifdef _GRAY_
  Mat img = imread (filename, CV_LOAD_IMAGE_GRAYSCALE);
  if (img.empty ())
    throw ERR_1;
  
  Mat src;
  img.convertTo (src, CV_64F);      
#else
  Mat img = imread(filename, CV_LOAD_IMAGE_COLOR);
  if (img.empty ())
    throw ERR_1;
  
  vector < cv::Mat> img_v;
  split(img, img_v);
  Mat img_g = img_v[1];
  Mat src;
  img_g.convertTo (src, CV_64F);
#endif // _GRAY_
#endif // _TEST_
}


void BlindDeblur::exe(void)
{
  exe(deblur->getImage());
}

void BlindDeblur::exe(const Mat img)
{
  //  param_t param;
  //  init_param (&param);

  /* for counter */
  calcTotalWorkload(&param);
  changeWeight(0);

#ifdef _TIME_
  START_TIME();
#endif

  try {
    /*-----------------------------------
     * Load blurred image
     *-----------------------------------*/

    // TODO :: COLOR 2 GRAY
    vector < cv::Mat> img_v;
    split(img, img_v);
    Mat img_g = img_v[1];
    Mat src;
    img_g.convertTo (src, CV_64F);

    /*
     * OpenCV to Eigen
     */
    MatrixXd Src (src.rows, src.cols);
    for (int x = 0; x < src.rows; x++)
      for (int y = 0; y < src.cols; y++)
	Src (x, y) = src.at <double>(x, y);
    
    /*
     * RGB
     */
    double _min, _max;
    normalizeImage (Src, _min, _max);
    
    /*
     * Get ROI
     */
    MatrixXd ROI = getROI (Src, param);
    
    if (param.MSlevels > MAX_LEVEL) {
      printf("ERROR: param.MSlevels (%d) is invalid. Less than %d\n",
	     param.MSlevels, MAX_LEVEL);
      throw ERR_2;
    }
    
    /*-----------------------------------
     * Create initial psf
     *-----------------------------------*/
    int h_rows = param.psf_size_x;
    int h_cols = param.psf_size_y;
    
    h_rows = h_rows >>(param.MSlevels -1);
    h_cols = h_cols >>(param.MSlevels -1);
    
    MatrixXd h = MatrixXd::Zero(h_rows, h_cols);  
    int cen_r = floor ((h.rows () + 1) / 2 - 1);
    int cen_c = floor ((h.cols () + 1) / 2 - 1);
    h (cen_r, cen_c) = 1.0;
    
    /*-----------------------------------
     * Estimate
     *-----------------------------------*/
    for (int L = 1; L <= param.MSlevels; L++) {
      if (checkInterrupt())
	throw INTR_1;

      if (L == param.MSlevels)
	changeWeight(16);
      else if (L == param.MSlevels -1)
	changeWeight(2);
      
#ifdef _TIME_
      SHOW_TIME("L=", L);
#endif
      /* Normalize psf     */
      matNormalize(h);
      
      /* Create ROI     */
      MatrixXd tmp = ROI;
      createROI(tmp, L, param);

#ifdef _VERBOSE0
      printf("createROI-out:\n");
#endif      
      /* Estimate PSF     */
      MatrixXcd cROI = MatrixXcd::Zero(tmp.rows(), tmp.cols());
      copy_mat_2_cmat(tmp, cROI, tmp.rows(), tmp.cols());
      PSFestimaLnoRgrad (h, cROI, param, L);
      
      /* double psf     */
      if (L != param.MSlevels)
	doublePSF(h);

    }
    
    /* Normalize PSF   */
    matNormalize(h);

    // centerize
    double g_x = 0.0; // center of gravity
    double g_y = 0.0;
    for (int x = 0; x < h.cols(); x++)
      for (int y = 0; y < h.rows(); y++) {
	g_x += x * h(x,y); g_y += y * h(x,y);
      }

    int shift_x = h.rows()/2 - (int)floor(g_x);
    int shift_y = h.cols()/2 - (int)floor(g_y);

    MatrixXd tmp = create_mat(h.rows(),h.cols());
    for (int x = 0; x < h.rows(); x++) 
      for (int y = 0; y < h.cols(); y++) {
	if (0 <= x+shift_x && x+shift_x < h.rows()
	    && 0 <= y+shift_y && y+shift_y < h.cols()) {
	  tmp(x+shift_x, y+shift_y) = h(x,y);
	}
	h(x,y) = 0.0;
      }
    
    for (int x = 0; x < h.rows(); x++) 
      for (int y = 0; y < h.cols(); y++)
	h(x,y) = tmp(x,y);
    
    /* Show psf   */
    show_psf(h);
    
    /*-----------------------------------
     * Deblurring
     *-----------------------------------*/

    // for counter
    changeWeight(23);

    double r_min, r_max, g_min, g_max, b_min, b_max;
    MatrixXd r_src, g_src, b_src;
    divide_rgb_images(img, r_src, g_src, b_src, r_min, r_max, g_min, g_max, b_min, b_max);
    
    MatrixXcd rU = create_complex_mat (r_src.rows (), r_src.cols ());
    MatrixXcd gU = create_complex_mat (g_src.rows (), g_src.cols ());
    MatrixXcd bU = create_complex_mat (b_src.rows (), b_src.cols ());

#ifdef _TIME_
    SHOW_TIME("fft", 1);
#endif
    fftCGSRaL (r_src, h, rU, param);
#ifdef _TIME_
    SHOW_TIME("RED", 1);
#endif
    fftCGSRaL (g_src, h, gU, param);
    
#ifdef _TIME_
    SHOW_TIME("GREEN", 1);
#endif
    fftCGSRaL (b_src, h, bU, param);
    
#ifdef _TIME_
    SHOW_TIME("RED", 1);
#endif

    /* Merge deblurred rgb images
     */
    Mat result;
    make_image(result, rU, gU, bU, r_min, r_max, g_min, g_max, b_min, b_max);
    
    /* Show deblur image
     */
    normalize (result, result, 0, 1, CV_MINMAX);

    /* Set Deblured Image and kernelImage
     */
    deblur->setDebluredImage(result);
    {
      Mat krnl(h.rows(), h.cols(), CV_64F);
      for (int x = 0; x < h.rows(); x++)
	for (int y = 0; y < h.cols(); y++)
	  krnl.at<double>(x,y) = h(x,y);
      normalize (krnl, krnl, 0, 1, CV_MINMAX);
      deblur->setKernelImage(krnl);
    }
  }
  catch (int err) {
    switch (err) {
    case ERR_1:
    ;
    break;
    case INTR_1:
    ;
    break;
    default:
      ;
    }
    //    throw err;
  }

}
  
/*-------------------------------------------------------
 * math functions
 *
 -------------------------------------------------------*/
double power (double x, double y)
{
  assert (x != 0.0);
  return exp (y * log (x));
}

double snx (setLnorm_t s)
{
  double y;
  double x = s.x;
  double alpha = s.alpha;
  double beta = s.beta;
  double q = s.q;
  
  y = -x + (1 - q) * q * alpha / beta * power (x, q - 1);
  return y;
}

double dsnx2 (setLnorm_t s)
{
  return snx (s);
}

double snx2 (setLnorm_t s)
{
  double y;
  double x = s.x;
  double alpha = s.alpha;
  double beta = s.beta;
  double q = s.q;

  y = -x * x / 2 + (1 - q) * alpha / beta * power (x, q);
  return y;
}

/*
 * Newton method
 *      fn : f(x)
 *      dfn : differential function of f(x)
 *      x0 : initial value
 *      eps1 : terminal condition 1（｜x(k+1)-x(k)｜＜eps1）
 *      eps2 : terminal condition 2（｜f(x(k))｜＜eps2）
 *      max : max trial number
 *      ind : actual trial number
 *      return : result
 */
double newton (double (*f) (setLnorm_t), double (*df) (setLnorm_t),
			    double alpha, double beta, double q,
			    double x0, double eps1, double eps2, int max, int *ind)
{
  double g, dg, x, x1;
  int sw;
  setLnorm_t sL;
  
  x1 = x0;
  x = x1;
  *ind = 0;
  sw = 0;

  sL.alpha = alpha;
  sL.beta = beta;
  sL.q = q;

  while (sw == 0 && *ind >= 0) {
      sw = 1;
      *ind += 1;
      sL.x = x1;
      g = (*f) (sL);

      if (fabs (g) > eps2) {
	if (*ind <= max) {
	      sL.x = x1;
	      dg = (*df) (sL);
	      if (fabs (dg) > eps2) {
		  x = x1 - g / dg;
		  if (fabs (x - x1) > eps1 && fabs (x - x1) > eps1 * fabs (x)) {
		      x1 = x;
		      sw = 0;
		  }
	      }
	      else
		*ind = -1;
	}
	else
	  *ind = -1;
      }
  }
  return x;
}

int BlindDeblur::sfix(const double d)
{
  if (d > 0)
    return floor(d);
  return ceil(d);
}


/*-------------------------------------------------------
 * image utils
 *
 *
 -------------------------------------------------------*/
double BlindDeblur::sumabs2 (const MatrixXcd & m)
{
  MatrixXcd tmp;
  tmp.array() = m.array().conjugate() * m.array();
  return real(tmp.sum());
}

void BlindDeblur::imresize (const MatrixXd &in, MatrixXd &out, const int sx, const int sy)
{
  int rsize = (int)in.rows();
  int csize = (int)in.cols();

  int type = INTER_LINEAR;
  //    int type = INTER_CUBIC; // INTER_NEAREST, INTER_LINER, INTER_AREA

  if (8 <= rsize && 8 <= csize)
    type = INTER_LANCZOS4;
  else if (4 <= rsize  && 4 <= csize)
    type = INTER_CUBIC;
  
  Mat src = Mat(rsize, csize, CV_64F);
  Mat dst = Mat(sx, sy, CV_64F);
  
  for (int x = 0; x < rsize; x++)
    for (int y = 0; y < csize; y++)
      src.at<double>(x, y) = in(x, y);
  
  resize(src, dst, Size(sx,sy), 0, 0, type);
  
  for (int x = 0; x < sx; x++)
    for (int y = 0; y < sy; y++)
      out(x, y) = dst.at<double>(x, y); 
}

void BlindDeblur::imresize (const MatrixXd &in, MatrixXd &out, const float m)
{
  int rsize = (int)in.rows();
  int csize = (int)in.cols();

  int type = INTER_LINEAR;
  //    int type = INTER_CUBIC; // INTER_NEAREST, INTER_LINER, INTER_AREA

  if (8 <= rsize && 8 <= csize)
    type = INTER_LANCZOS4;
  else if (4 <= rsize  && 4 <= csize)
    type = INTER_CUBIC;
  
  Mat src = Mat(rsize, csize, CV_64F);
  Mat dst = Mat(sfix(rsize*m), sfix(csize*m), CV_64F);
  
  for (int x = 0; x < rsize; x++)
    for (int y = 0; y < csize; y++)
      src.at<double>(x, y) = in(x, y);
  
  resize(src, dst, Size(), m, m, type);
  
  for (int x = 0; x < sfix(rsize*m); x++)
    for (int y = 0; y < sfix(csize*m); y++)
      out(x, y) = dst.at<double>(x, y); 
}

void BlindDeblur::imresize (const MatrixXd &in, MatrixXd &out)
{
  //  imresize (in, out, (float)2.0);
  imresize (in, out, 2*in.rows(), 2*in.cols());
}

void BlindDeblur::normalizeImage (MatrixXd & m, double &min, double &max)
{
  min = m.minCoeff ();
  max = m.maxCoeff ();

  double scale = max - min;
  if (scale < 0.0)
    scale *= -1.0;

  assert ((scale) != 0.0);
  m.array () -= min;
  m.array () /= scale;
}

/*-------------------------------------------------------
 * asetLnorm
 *
 -------------------------------------------------------*/
void BlindDeblur::aLn (MatrixXd & V, const MatrixXd & DU, const MatrixXd & normDU, const vu_star_t r)
{
  double u_star = r.u_star;
  double v_star = r.v_star;
  double k = u_star - v_star;

  assert(u_star > 0);
  MatrixXd mask = MatrixXd::Zero (V.rows(), V.cols());
  mask = (normDU.array() > u_star).select(1.0, normDU*0.0);
  V.array() = DU.array() * (normDU.array() -k) * normDU.array().inverse() * mask.array();
}

void BlindDeblur::asetLnorm (const double q, double alpha, double beta, vu_star_t * r)
{
  if (1.0 == q) {
      r->v_star = 0.0;
      r->u_star = (double)(alpha / beta);
  }
  else if (q == 0.0) {
    r->v_star = sqrt (2 * alpha / beta);
    r->u_star = sqrt (2 * alpha / beta);
  }
  else {
    double eps1, eps2, x0;
    int max, ind;
    eps1 = 2.0e-16;
    eps2 = 2.0e-16;
    max = 100;
    x0 = 0.1;
    //    double leftmarker;
    //    leftmarker = newton(snx, dsnx, alpha, beta, q, x0, eps1, eps2, max, &ind);
    //    double tmp = newton(snx2, dsnx2, alpha, beta, q, x0, leftmarker, leftmarker, max, &ind);
    r->v_star =
      newton (snx2, dsnx2, alpha, beta, q, x0, eps1, eps2, max, &ind);
    r->u_star = r->v_star + alpha / beta * q * power (r->v_star, q - 1);
  }
}


/*-------------------------------------------------------
 * PSFestimaLnoRgrad
 *
 *
 -------------------------------------------------------*/

void BlindDeblur::get_margins_mat (const MatrixXd & m, int &m_left, int &m_right, int &m_top,
				   int &m_bottom,
				   const double threshold)
{
  // for H
  m_left = 1;
  m_right = m.rows();
  m_top = 1;
  m_bottom = m.cols();

  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)
      if (m (x, y) > threshold) {
	m_top = x+1;
	goto out1;
      }
 out1:
  for (int x = m.rows () - 1; 0 <= x; x--)
    for (int y = 0; y < m.cols (); y++)
      if (m (x, y) > threshold) {
	m_bottom = x + 1;
	goto out2;
      }
 out2:
  for (int y = 0; y < m.cols (); y++)
    for (int x = 0; x < m.rows (); x++)
      if (m (x, y) > threshold) {
	m_left = y+1;
	goto out3;
      }
 out3:
  for (int y = m.cols () - 1; 0 <= y; y--)
    for (int x = 0; x < m.rows (); x++)
      if (m (x, y) > threshold) {
	m_right = y+1;
	goto out4;
      }
 out4:
  ;
}

void BlindDeblur::make_mask_mat (const MatrixXd & m, MatrixXd & mask, const double threshold)
{
  // for H
  mask = (m.array() > threshold).select(1.0, m*0.0);
}

void BlindDeblur:: mat2gray (MatrixXd & m)
{
  // for H
  double min, max;
  normalizeImage (m, min, max);
}

void BlindDeblur::bwmorph_clean (MatrixXd & mask)
{
  int m_x = (int)mask.rows();
  int m_y = (int)mask.cols();

  // for H
  MatrixXd temp = MatrixXd::Zero (m_x + 2, m_y + 2);
  for (int x = 0; x < m_x; x++)
    for (int y = 0; y < m_y; y++)
      temp (x + 1, y + 1) = mask (x, y);

  for (int x = 1; x <= m_x; x++)
    for (int y = 1; y <= m_y; y++)
      if (
	  temp (x - 1, y - 1) == 0.0  && 
	  temp (x - 1, y) == 0.0 &&
	  temp (x - 1, y + 1) == 0.0 && 
	  temp (x,     y - 1) == 0.0 &&
	  temp (x,     y + 1) == 0.0 &&
	  temp (x + 1, y - 1) == 0.0 &&
	  temp (x + 1, y) == 0.0 &&
	  temp (x + 1, y + 1) == 0.0)
       	mask (x - 1, y - 1) = 0.0;
}


/* 
 * Vh(Vh<0)=0; 
 * Vh(hsize(1)+1:end,:,:) = 0; Vh(1:hsize(1),hsize(2)+1:end,:) = 0;
 */
void BlindDeblur::set_Vh (MatrixXd & vh, const int rows, const int cols)
{
#ifdef _TEST
  print_mat2(vh, 0, "set_Vh", "Vh-init");
#endif

  MatrixXd tmp = MatrixXd::Zero (rows, cols);
  tmp = vh.block(0,0,rows-1, cols-1);
  tmp = (tmp.array() < 0.0).select(0.0, tmp);

  vh.setZero();
  vh.block(0,0,rows-1,cols-1).array() = tmp.array();

#ifdef _TEST
  print_mat2(vh, 0, "set_Vh", "Vh-last");
#endif
}

/*
 * centerPSF
 */
void BlindDeblur::centerPSF (MatrixXd & H, const param_t param)
{
  int m_left, m_right, m_top, m_bottom;
  
  mat2gray (H);
  MatrixXd mask = MatrixXd::Zero (H.rows (), H.cols ());
  make_mask_mat (H, mask, (double) param.centering_threshold);
  bwmorph_clean (mask);
  get_margins_mat (mask, m_left, m_right, m_top, m_bottom, (double) param.centering_threshold);

  int topleft_x = m_top;
  int topleft_y = m_left;

  int begin_x = smax(topleft_x, 1) -1;
  int end_x = smin(topleft_x + H.rows() -1, H.rows()) -1;
  int begin_y = smax(topleft_y, 1) -1;
  int end_y = smin(topleft_y + H.cols() -1, H.cols()) -1;

  int shift_x = (begin_x + end_x - H.rows())/2;
  int shift_y = (begin_y + end_y - H.cols())/2;

#ifdef _DEBUG01
  printf("CENTERPSF: (%d - %d, %d - %d) shift(%d,%d):%d,%d\n", begin_x, end_x, begin_y, end_y, shift_x, shift_y, H.rows(), H.cols());
  fflush(stdout);
#endif

#ifdef _DEBUG01
  printf("top=%d bottom=%d left=%d right=%d\n", m_top,m_bottom, m_left, m_right);
  printf("topleft(%d,%d) begin(%d,%d) end(%d,%d)\n", topleft_x, topleft_y,
	 begin_x, begin_y, end_x, end_y);
#endif

  MatrixXd tmp = MatrixXd::Zero (H.rows (), H.cols ());

  /*
  for (int x = begin_x - shift_x; x <= end_x - shift_x; x++)
    for (int y = begin_y - shift_y; y <= end_y - shift_y; y++)
      if (H(x,y) <= param.centering_threshold)
	tmp(x,y) = 0.0;
	else
	tmp (x, y) = H (x-shift_x, y-shift_y);
  */
  
  for (int x = begin_x; x <= end_x; x++)
    for (int y = begin_y; y <= end_y; y++)
      if (H(x,y) <= param.centering_threshold)
	tmp(x,y) = 0.0;
      else
	tmp (x, y) = H (x, y);
  // 	tmp (x-shift_x, y-shift_y) = H (x, y);
  
  tmp.swap (H);
  double s = H.sum();
  if (s != 0.0)
    H.array() /= s;

}


void BlindDeblur::PSFestimaLnoRgrad (MatrixXd & h, MatrixXcd & ROI, param_t param, const int L)
{
  try {
    
    double gamma = param.gamma;
    int hsize_r = h.rows ();
    int hsize_c = h.cols ();
    int gsize_r = ROI.rows ();
    int gsize_c = ROI.cols ();
    int usize_r = gsize_r;
    int usize_c = gsize_c;
    
    /*
     * Matrix
     */
    
    /* common
     */
    MatrixXd U = create_mat(gsize_r, gsize_c);
    MatrixXcd FU = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FUx = create_complex_mat (usize_r, usize_c);
    MatrixXcd FUy = create_complex_mat (usize_r, usize_c);
    MatrixXcd FDx = create_complex_mat (usize_r, usize_c);
    MatrixXcd FDy = create_complex_mat (usize_r, usize_c);
    MatrixXcd DTD = create_complex_mat (usize_r, usize_c);
    MatrixXcd b = create_complex_mat (usize_r, usize_c);
    MatrixXcd FeGu = create_complex_mat (usize_r, usize_c);
    MatrixXcd FeGx = create_complex_mat (usize_r, usize_c);
    MatrixXcd FeGy = create_complex_mat (usize_r, usize_c);
    MatrixXd Vx = create_mat (usize_r, usize_c);
    MatrixXd Vy = create_mat (usize_r, usize_c);
    MatrixXd Vh = create_mat (usize_r, usize_c);
    MatrixXd Bx = create_mat (usize_r, usize_c);
    MatrixXd By = create_mat (usize_r, usize_c);
    MatrixXd Bh = create_mat (usize_r, usize_c);
    
    /* for Ustep
     */
    MatrixXcd FHS = create_complex_mat (usize_r, usize_c);
    MatrixXcd FHTH = create_complex_mat (usize_r, usize_c);
    MatrixXcd FGs = create_complex_mat (usize_r, usize_c);
    MatrixXcd FUp = create_complex_mat (usize_r, usize_c);
    
    MatrixXd xD = create_mat (usize_r, usize_c);
    MatrixXd yD = create_mat (usize_r, usize_c);
    MatrixXd xDm = create_mat (usize_r, usize_c);
    MatrixXd yDm = create_mat (usize_r, usize_c);
    MatrixXd nDm = create_mat (usize_r, usize_c);
    
    /* for Hstep
     */
    MatrixXcd FUD = create_complex_mat (usize_r, usize_c);
    MatrixXcd FUTU = create_complex_mat (usize_r, usize_c);
    MatrixXcd FH = create_complex_mat (usize_r, usize_c);
    MatrixXcd FHp = create_complex_mat (usize_r, usize_c);
    MatrixXd hI = create_mat (usize_r, usize_c);
    MatrixXd hIm = create_mat (usize_r, usize_c);
    MatrixXd nIm = create_mat (usize_r, usize_c);
    
    /* utils
     */
    MatrixXcd tmp1 = create_complex_mat (usize_r, usize_c);
    MatrixXcd tmp2 = create_complex_mat (usize_r, usize_c);
    MatrixXd tmpr1 = create_mat (usize_r, usize_c);
    MatrixXd tmpr2 = create_mat (usize_r, usize_c);
    

    if (checkInterrupt())
      throw INTR_1;

    /*
     *
     */
    dft2d (U, FU);
    
    complex < double >c;
    c = complex < double >(1.0, 0.0);
    FDx(0, 0) = c;
    FDy(0, 0) = c;
    
    c = complex < double >(-1.0, 0.0);
    FDx(0, 1) = c;
    FDy(1, 0) = c;
    
    dft2d (FDx, FDx);
    dft2d (FDy, FDy);
    
    /* DTD = conj(FDx).*FDx + conj(FDy).*FDy */
    DTD.array() = FDx.array().conjugate() * FDx.array() + FDy.array().conjugate() * FDy.array();

    /*
     * edgetaper: ROI
     */
#ifdef _EDGETAPER
    MatrixXcd eG = create_complex_mat (usize_r, usize_c);
    edgetaper(eG, ROI);
    dft2d (eG, FeGu);
#else
    dft2d (ROI, FeGu);
#endif
    
    FeGx.array() = FDx.array() * FeGu.array();  //  FeGx = FDx .* FeGu;
    FeGy.array() = FDy.array() * FeGu.array();  //  FeGy = FDy .* FeGu;
    
    /*
     *
     */
    for (int ml = 1; ml <= param.maxiter; ml++)
      {
	if (checkInterrupt())
	  throw INTR_1;


#ifdef _TIME_
	SHOW_TIME("PSF=", ml);
#endif


#ifdef _VERBOSE0
	printf("PSFestimaLnoRgrad:L=%d loop=%d\n", L, ml); fflush(stdout);
#endif
	
	Ustep (param.Lp, h, U, FeGu, FU, FUx, FUy, FDx, FDy, Vx, Vy, Bx, By,
	       DTD, FHS, FHTH, FGs, FUp, b, xD, yD, xDm, yDm, nDm, tmp1, tmp2,
	       tmpr1, tmpr2, usize_r, usize_c, param, gamma, L, ml);
	
	Hstep (1.0, h, FeGx, FeGy, FUx, FUy, Vh, Bh,
	       FHp, hI, hIm, nIm, FUD, FUTU, FH, b, tmp1, tmp2, tmpr1,
	       usize_r, usize_c, hsize_r, hsize_c, param, gamma, L, ml);
	gamma *= 1.5;
      }
#ifdef _VERBOSE0
    printf("\n");
#endif
    
    centerPSF (h, param);
    
#ifdef _VERBOSE1
    std::cout << "HI =" << endl << h << endl;
#endif
    
  }
  catch (int err) {
    throw err;
  }
}


/*
 * Ustep
 */
void BlindDeblur::Ustep (const double q,
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
			 const int usize_r, const int usize_c, param_t param, double gamma, const int L, const int ml)
{
  double beta = param.beta_u;
  double alpha = param.alpha_u;
  double ccreltol = param.ccreltol;

  try {

    /* dft(FHS, usize_h, usize_w); */
    copy_mat_2_mat_zeros (H, tmpr1, H.rows (), H.cols ());
    dft2d (tmpr1, FHS);
    
    /* FHTH = conj(FHS).* FHS */
    FHTH.array() = FHS.array().conjugate() * FHS.array();
    /* FGs = sum(conj(FHS).* FeGu, 3) */  
    FGs.array() = FHS.array().conjugate() * FeGu.array();
    
#ifdef _VERBOSE0
    printf("USTEP loop:");
#endif
    
    for (int i = 1; i <= param.maxiter_u; i++)
      {
	if (checkInterrupt())
	  throw INTR_1;	

	/* for counter */
	countUp();	

#ifdef _VERBOSE0
	printf("%3d ", i); fflush(stdout);
#endif
	
	FUp = FU;      //      copy_cmat (FU, FUp, usize_r, usize_c);
	
	/* b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By)); */
	tmpr1.array() = Vx.array() + Bx.array();
	
	dft2d (tmpr1, tmp1);
	
	tmpr1.array() = Vy.array() + By.array();
	
	dft2d (tmpr1, tmp2);
	
	b.array() = FDx.array().conjugate() * tmp1.array() + FDy.array().conjugate() * tmp2.array();
	
	//      b = (beta / gamma)*b + FGs;
	b.array() = beta/gamma*b.array() + FGs.array();
	
	/* FU = b ./ (FHTH + beta/gamma*DTD); */
	tmp1.array() = beta/gamma*DTD.array() + FHTH.array();
	
	FU.array() = b.array() * tmp1.array().inverse(); 
	
	FUx.array() = FDx.array() * FU.array();      /* FUx = FDx .* FU */
	FUy.array() = FDy.array() * FU.array();      /* FUy = FDy .* FU */
	
	/* xD = real(iff2(FUx)) */
	/* yD = real(iff2(FUy)) */
	idft2d (FUx, xD);
	idft2d (FUy, yD);
	
	xDm.array() = xD.array() - Bx.array();
	yDm.array() = yD.array() - By.array();
	
	/* nDm = sqrt(xDm.^2 + yDm.^2) */
	nDm.array() = xDm.array() * xDm.array() + yDm.array() * yDm.array();
	nDm.array() = nDm.array().sqrt();
	
	/* Vx = Pr.fh(xDm, nDm) */
	/* Vy = Pr.fh(yDm, nDm) */
	vu_star_t vu;
	asetLnorm (q, alpha, beta, &vu);
	aLn (Vx, xDm, nDm, vu);
	aLn (Vy, yDm, nDm, vu);
	
	Bx.array() = Bx.array() + (Vx.array() - xD.array());
	By.array() = By.array() + (Vy.array() - yD.array());
	
	double relcon = sqrt (sumabs2 (FUp.array() - FU.array())) / sqrt (sumabs2 (FU));
	
	if (relcon < ccreltol)
	  break;
      }
#ifdef _VERBOSE0
    printf("\n");
#endif
    
    
    /* U = real(ifft2(FU)) */
    idft2d (FU, U);
    
  }
  catch(int err) {
    throw err;
  }
}

/*
 * Hstep
 */
void BlindDeblur::Hstep (double q,
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
			 param_t param, double gamma, const int L, const int ml)
{
  double relcon;
  double ccreltol = param.ccreltol;
  double beta = param.beta_h;
  double alpha = param.alpha_h;

  try {

    /* FUD = FeGx.*conj(FUx) + FeGy.*conj(FUy); */
    FUD.array() = FeGx.array() * FUx.array().conjugate() + FeGy.array() * FUy.array().conjugate();
    /* FUTU = conj(FUx).*FUx + conj(FUy).*FUy; */
    FUTU.array() = FUx.array().conjugate() * FUx.array()  + FUy.array().conjugate() * FUy.array();
    
    /* FH = fft2(H,usize(1),usize(2)); */
    copy_mat_2_cmat_zeros (H, tmp1, H.rows (), H.cols ());
    dft2d (tmp1, FH);
    
    /*   tmp1 = (FUTU + beta/gamma) */
    tmp1.real().array() = FUTU.real().array() + (beta / gamma);
    
#ifdef _VERBOSE0
    printf ("HSTEP: loop=");
#endif
    
    for (int i = 1; i <= param.maxiter_h; i++)
      {
	if (checkInterrupt())
	  throw INTR_1;

	/* for counter */
	countUp();
	
#ifdef _VERBOSE0
	printf ("%3d", i); fflush(stdout);
#endif
	
	FHp = FH;
	
	/* b = beta/gamma*fft2(Vh+Bh) + FUD */
	tmpr1.array() = Vh.array() + Bh.array();
	
	dft2d (tmpr1, b);

	//b = (beta / gamma)* b + FUD;
	b.array() = (beta / gamma)* b.array() + FUD.array();
	
	/* FH = b./(FUTU + beta/gamma) = b./tmp1 */
	FH.array() = b.array() * tmp1.array().inverse();
	
	/* relcon = sqrtsum(abs(FHp(:) - FH(:)).^2) / sqrt(sum(abs(FH(:)).^2)) */
	relcon = sqrt (sumabs2 (FHp.array() - FH.array())) / sqrt (sumabs2 (FH));
	idft2d (FH, hI);
	
	hIm.array() = hI.array() - Bh.array();
	
	/* nIm = abs(hIm); */
	nIm.array() = hIm.array().abs();
	
	/* Vh = Pr.fh(hIm, nIm); */
	vu_star_t vu;
	asetLnorm (q, alpha, beta, &vu);
	aLn (Vh, hIm, nIm, vu);
	/* 
	 * Vh(Vh<0);
	 * Vh(hsize(1)+1:end,:,:) = 0; 
	 * Vh(1:hsize(1),hsize(2)+1:end,:) = 0;
	 */
	set_Vh (Vh, hsize_r, hsize_c);
	Bh.array() = Bh.array() + (Vh.array() - hI.array());
	
	/*  H = hI(1:hsize(1),1:hsize(2),:); */
	H = hI.block(0, 0, hsize_r, hsize_c);
	
	if (relcon < ccreltol)
	  break;
      }
    
#ifdef _VERBOSE0
    printf("\n");
#endif
    
  }
  catch(int err) {
    throw err;
  }
}

/*-------------------------------------------------------
 * fftCGSRaL
 *
 -------------------------------------------------------*/
void BlindDeblur::set_vrange (const MatrixXd & m, vrange_t & vrange)
{
  vrange.min = m.minCoeff ();
  vrange.max = m.maxCoeff ();
}

void BlindDeblur::uConstr (MatrixXcd & m, const vrange_t vrange)
{
  for (int x = 0; x < m.rows (); x++)
    for (int y = 0; y < m.cols (); y++)  {
      if (real (m (x, y)) < vrange.min)
	m (x, y) = complex < double >(vrange.min, 0.0);
      if (real (m (x, y)) > vrange.max)
	m (x, y) = complex < double >(vrange.max, 0.0);
    }
}

void BlindDeblur::fftCGSRaL (MatrixXd & G, MatrixXd & H, MatrixXcd & U, const param_t param)
{
  try {
    int maxiter = param.maxiter_u;
    double alpha = param.alpha_u;
    double ccreltol = param.ccreltol;
    double gamma = param.gamma_nonblind;
    double beta = param.beta_u_nonblind;
    double Lp = param.Lp_nonblind;
    
    double gsize_r = G.rows ();
    double gsize_c = G.cols ();
    double hsize_r = H.rows ();
    double hsize_c = H.cols ();
    
    /*
     */
    MatrixXcd FU = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FDx = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FDy = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FH = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FHTH = create_complex_mat (gsize_r, gsize_c);
    
    MatrixXcd tmp1 = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd tmp2 = create_complex_mat (gsize_r, gsize_c);
    MatrixXd tmpr1 = create_mat (gsize_r, gsize_c);
    MatrixXd tmpr2 = create_mat (gsize_r, gsize_c);
    
    vrange_t vrange;
    set_vrange (G, vrange);
    
    MatrixXcd hshift = create_complex_mat (hsize_r, hsize_c);
    {
      int cen_r, cen_c;
      cen_r = floor ((H.rows () + 1) / 2);
      cen_c = floor ((H.cols () + 1) / 2);
      hshift (cen_r, cen_c) = complex < double >(1.0, 0.0);
    }
    
    complex < double >c;
    c = complex < double >(1.0, 0.0);
    FDx(0, 0) = c;
    FDy(0, 0) = c;
    
    c = complex < double >(-1.0, 0.0);
    FDx(0, 1) = c;
    FDy(1, 0) = c;
    
    dft2d (FDx, FDx);
    dft2d (FDy, FDy);
    
    /* FH = conj(fft2(hshift,usize)) .* fft2(H, usize) */
    copy_cmat_2_cmat_zeros (hshift, tmp1, hshift.rows (), hshift.cols ());
    copy_mat_2_cmat_zeros (H, tmp2, hshift.rows (), hshift.cols ());
    dft2d (tmp1, tmp1);
    dft2d (tmp2, tmp2);
    FH.array() = tmp1.array().conjugate() * tmp2.array();
    FHTH.array() = FH.array().conjugate() * FH.array();
    
    MatrixXcd FGu = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd FGs = create_complex_mat (gsize_r, gsize_c);
    
    MatrixXcd cG = create_complex_mat (gsize_r, gsize_c);
    copy_mat_2_cmat (G, cG, G.rows (), G.cols ());
    /*
     * edgetaper
     */
#ifdef _EDGETAPER
    MatrixXcd eG = create_complex_mat (gsize_r, gsize_c);
    edgetaper(eG, cG);
    dft2d (eG, FGu);
#else
    dft2d (cG, FGu);
#endif
    
    FGs.array() = FH.array().conjugate() * FGu.array();
    
    MatrixXcd DTD = create_complex_mat (gsize_r, gsize_c);
    DTD.array() = FDx.array().conjugate() * FDx.array() + FDy.array().conjugate() * FDy.array();
    
    /*
     *
     */
    MatrixXcd FUp = create_complex_mat (gsize_r, gsize_c);
    MatrixXcd b = create_complex_mat (gsize_r, gsize_c);
    
    MatrixXd xD = create_mat (gsize_r, gsize_c);
    MatrixXd yD = create_mat (gsize_r, gsize_c);
    MatrixXd xDm = create_mat (gsize_r, gsize_c);
    MatrixXd yDm = create_mat (gsize_r, gsize_c);
    MatrixXd nDm = create_mat (gsize_r, gsize_c);
    
    MatrixXd Bx = create_mat (gsize_r, gsize_c);
    MatrixXd By = create_mat (gsize_r, gsize_c);
    MatrixXd Vx = create_mat (gsize_r, gsize_c);
    MatrixXd Vy = create_mat (gsize_r, gsize_c);
    
#ifdef _VERBOSE0
    printf ("FFT-LOOP=");
#endif
    
    for (int i = 1; i <= maxiter; i++)
      {
	if (checkInterrupt())
	  throw INTR_1;	

	/* for counter */
	countUp();
	
#ifdef _VERBOSE0
	printf ("%3d", i);
#endif
	
	/* FUp = FU; */
	FUp = FU;
	
	/*  b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By)); */
	tmpr1.array() = Vx.array() + Bx.array();      dft2d (tmpr1, tmp1);
	tmpr1.array() = Vy.array() + By.array();      dft2d (tmpr1, tmp2);
	b.array() = FDx.array().conjugate() * tmp1.array() + FDy.array().conjugate() * tmp2.array();
	
	//      b = b * beta / gamma + FGs;
	b.array() = beta/gamma*b.array() + FGs.array();
	
	/* FU = b./(FHTH + beta/gamma*DTD); */
	tmp1.array() = beta/gamma*DTD.array() + FHTH.array();
	
	FU.array() = b.array() * tmp1.array().inverse();
	
	/* xD = real(ifft2(FDx.*FU)); */
	/* yD = real(ifft2(FDy.*FU)); */
	tmp1.array() = FDx.array() * FU.array();
	idft2d (tmp1, xD);
	tmp2.array() = FDy.array() * FU.array();
	idft2d (tmp2, yD);
	xDm.array() = xD.array() - Bx.array();
	yDm.array() = yD.array() - By.array();
	
	/* nDm = repmat(sqrt(sum(xDm.^2,3) + sum(yDm.^2,3)),[1 1 usize(3)]); */
	nDm.array() = xDm.array() * xDm.array() + yDm.array() * yDm.array();
	nDm.array() = nDm.array().sqrt();
	
	/* Pr = asetupLnormPrior(Lp,alpha,beta); */
	vu_star_t vu;
	asetLnorm (Lp, alpha, beta, &vu);
	
	/* Vy = Pr.fh(yDm,nDm); */
	/* Vx = Pr.fh(xDm,nDm); */
	aLn (Vx, xDm, nDm, vu);
	aLn (Vy, yDm, nDm, vu);
	
	/* Bx = Bx + Vx - xD; */
	/* By = By + Vy - yD; */
	Bx.array() = Bx.array() + (Vx.array() - xD.array());
	By.array() = By.array() + (Vy.array() - yD.array());
	
	/* relcon = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2)); */
	double relcon = sqrt (sumabs2 (FUp - FU)) / sqrt (sumabs2 (FU));
	
	if (relcon < ccreltol)
	  break;
      }
#ifdef _VERBOSE0
    printf("\n");
#endif
    
    /* U = real(ifft2(FU)); */
    idft2d (FU, U);
    
    uConstr (U, vrange);
    
  }
  catch(int err) {
    throw err;
  }
}

/*-------------------------------------------------------
 * utils
 *
 *
 -------------------------------------------------------*/

MatrixXd BlindDeblur::getROI (const MatrixXd & m, const param_t param)
{
  int m_r = m.rows ();
  int m_c = m.cols ();
  int margin_r = 0;
  int margin_c = 0;

  if (m_r > param.maxROIsize_r) {
    margin_r = floor ((m_r - param.maxROIsize_r) / 2);
    m_r = param.maxROIsize_r;
  }
  if (m_c > param.maxROIsize_c) {
    margin_c = floor ((m_c - param.maxROIsize_c) / 2);
    m_c = param.maxROIsize_c;
  }
  return m.block (margin_r, margin_c, m_r + margin_r, m_c + margin_c);
}


void BlindDeblur::createROI(MatrixXd &tmp, const int L, const param_t param)
{
  if (L != param.MSlevels) {
    int tmp_x = (int)tmp.rows();
    int tmp_y = (int)tmp.cols();
    float m = 1.0;
    for (int i = L; i < param.MSlevels; i++) {
      tmp_x *= 0.5;	tmp_y *= 0.5;
      m *= 0.5;
    }
    MatrixXd tmp2 = MatrixXd::Zero(tmp_x, tmp_y);
    imresize(tmp, tmp2, m);
    
#ifdef _DEBUG1
    printf("ROI(%d), %d, %d, %d %d\n", 
	   L, (int)tmp.rows(), (int)tmp.cols(), (int)tmp2.rows(), (int)tmp2.cols());
#endif

    tmp.resize(tmp_x, tmp_y);
    tmp2.swap(tmp);
  }
}
    
void BlindDeblur::doublePSF(MatrixXd & h)
{
  int hx = h.rows();
  int hy = h.cols();
  MatrixXd hi = MatrixXd::Zero(hx, hy);
  hi.array() = h.array();
  h.resize(hx*2, hy*2);
  //imresize(hi, h);
  imresize(hi, h, hx*2, hy*2);
}

/*-------------------------------------------------------
 * display utils
 *
 *
 -------------------------------------------------------*/
void BlindDeblur::divide_rgb_images(Mat color_img,
				    MatrixXd & r_src, MatrixXd & g_src, MatrixXd & b_src,
				    double & r_min, double & r_max,
				    double & g_min, double & g_max,
				    double & b_min, double & b_max)
{
  vector < cv::Mat> img_planes;
  split(color_img, img_planes);

  Mat r_img = img_planes[0];
  Mat g_img = img_planes[1];
  Mat b_img = img_planes[2];
  Mat r_img2, g_img2, b_img2;
  r_img.convertTo (r_img2, CV_64F);
  g_img.convertTo (g_img2, CV_64F);
  b_img.convertTo (b_img2, CV_64F);
  
  r_src = create_mat(r_img.rows, r_img.cols);
  g_src = create_mat(g_img.rows, g_img.cols);
  b_src = create_mat(b_img.rows, b_img.cols);
  
  for (int x = 0; x < r_img2.rows; x++)
    for (int y = 0; y < r_img2.cols; y++) {
      r_src (x, y) = r_img2.at<double>(x, y);
      g_src (x, y) = g_img2.at<double>(x, y);
      b_src (x, y) = b_img2.at<double>(x, y);
    }

  normalizeImage (r_src, r_min, r_max);
  normalizeImage (g_src, g_min, g_max);
  normalizeImage (b_src, b_min, b_max);

}


void BlindDeblur::make_image(Mat & result,
			     const MatrixXcd rU, const MatrixXcd gU, const MatrixXcd bU, 
			     const double r_min, const double r_max,
			     const double g_min, const double g_max,
			     const double b_min, const double b_max)
{
  Mat out[3] = Mat::zeros(rU.rows (), rU.cols (), CV_64F);
  for (int x = 0; x < rU.rows (); x++)
    for (int y = 0; y < rU.cols (); y++) {
      out[0].at < double >(x, y) = real (rU (x, y)) * (r_max - r_min) + r_min;
      out[1].at < double >(x, y) = real (gU (x, y)) * (g_max - g_min) + g_min;
      out[2].at < double >(x, y) = real (bU (x, y)) * (b_max - b_min) + b_min;
    }
  
  vector < cv::Mat > ret;
  ret.push_back(out[0]);
  ret.push_back(out[1]);
  ret.push_back(out[2]);
  cv::merge(ret, result);

}

void BlindDeblur::show_psf(const MatrixXd &h)
{
  double s = h.sum ();
  Mat kernel = Mat::zeros(h.rows(), h.cols(), CV_64F);
  for (int x = 0; x < h.rows(); x++)
    for (int y = 0; y < h.cols(); y++)
      kernel.at<double>(x,y) = h(x,y) / s;
  
  normalize (kernel, kernel, 0, 1, CV_MINMAX);
  
  /* Set kernel Image
   */
  if (deblur != NULL) {
    deblur->setKernelImage(kernel);
  }
}


void BlindDeblur::matNormalize(MatrixXd & h)
{
  double s;
  s = h.sum ();
  if (s > 0.0)
    h.array () /= s;
}


void BlindDeblur::sharp_filter(Mat &result, const param_t param)
{
  double f = param.filter;
  
  if (f != 0.0) {
    Mat kernel = Mat::zeros(3, 3, CV_64F);
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
	kernel.at<double>(x,y) = -f / 9.0;
    kernel.at<double>(1,1) = 1+(8*f)/ 9.0;
    
    filter2D(result,result, -1, kernel);
  }
}

/*-------------------------------------------------------
 * param
 *
 *
 -------------------------------------------------------*/

void BlindDeblur::init_param (param_t * param)
{
  param->mode = 0; // dummy

  param->psf_size_x = 256;  // 32
  param->psf_size_y = 256;  // 32

  param->maxROIsize_r = 1024;
  param->maxROIsize_c = 1024;

  param->MSlevels = 4;   // 4

  param->gamma = 1e2; // 1e2
  param->Lp = 0.4;    // 0.3　　　// if exceed 0.5,  K＝0

  param->beta_h = 1e4 * param->gamma;  // 1e4*gamma  
  param->alpha_h = 1e1 * param->gamma; // 1e1*gamma　
  param->centering_threshold = (double) 30/255; // 20/255

  param->beta_u = 1e0 * param->gamma;  // 1e0 * gamma  
  param->alpha_u = 1e-2 * param->gamma; // 1e-2 * gamma

  param->gamma_nonblind = 2e3 * param->gamma; // 2e1 * gamma
  param->beta_u_nonblind = 1.0 * param->gamma_nonblind; // 1e-2 * gamma
  param->Lp_nonblind = 0.0;  // 1.0

  param->maxiter_u = 7;	//10
  param->maxiter_h = 7;	//10
  param->maxiter = 7;		//10
  param->ccreltol = 1e-3;  // 1e-3

  param->filter = 3.5;
}

void BlindDeblur::init_param (void)
{
  param.psf_size_x = 32;
  param.psf_size_y = 32;

  param.maxROIsize_r = 1024;
  param.maxROIsize_c = 1024;

  param.MSlevels = 4;   // 4

  param.gamma = 1e2; // 1e2
  param.Lp = 0.3;    // 0.3　　　// if exceed 0.5, K＝0

  param.beta_h = 1e4 * param.gamma;  // 1e4*gamma  
  param.alpha_h = 1e1 * param.gamma; // 1e1*gamma　
  param.centering_threshold = (double) 30/255; // 20/255

  param.beta_u = 1e0 * param.gamma;  // 1e0 * gamma  
  param.alpha_u = 1e-2 * param.gamma; // 1e-2 * gamma　

  param.gamma_nonblind = 2e3 * param.gamma; // 2e1 * gamma
  param.beta_u_nonblind = 1.0 * param.gamma_nonblind; // 1e-2 * gamma
  param.Lp_nonblind = 0.0;  // 1.0

  param.maxiter_u = 10;	//10
  param.maxiter_h = 10;	//10
  param.maxiter = 5;		//10
  param.ccreltol = 1e-3;  // 1e-3

  param.filter = 3.5;
}

void BlindDeblur::setPsfSize(const int size)
{
  param.psf_size_x = size;
  param.psf_size_y = size;
}

void BlindDeblur::setMode(const int mode)
{
  param.mode = mode;

  switch (param.mode) {
  case DEBLUR_LOW:
    param.MSlevels = 1;
    break;
  case DEBLUR_MEDIUM:
    param.MSlevels = 2;
    break;
  default:
    param.MSlevels = 4;
  }

}

/*
 *
 */
void BlindDeblur::calcTotalWorkload(param_t * param) {
  count = 0;
  totalWorkload = 0;
  counter = 0;

  weight = 16;
  totalWorkload += weight * param->maxiter_u * param->maxiter;
  totalWorkload += weight * param->maxiter_h * param->maxiter;

  if (param->MSlevels > 1) {
    weight = 2;
    totalWorkload += weight * param->maxiter_u * param->maxiter;
    totalWorkload += weight * param->maxiter_h * param->maxiter;
  }
  weight = 23;
  totalWorkload += 3 * weight * param->maxiter_u; // 3 colors

  // reset weight value
  weight = 1;
};

void BlindDeblur::countUp(void) {
  counter += weight;
  count = floor((counter*100)/totalWorkload); // 0 to 100
};

void BlindDeblur::changeWeight(const int weight) {
  this->weight = weight;
};

// EOF
