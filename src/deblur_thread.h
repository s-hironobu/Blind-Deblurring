#ifndef _DEBLUR_THREAD_
#define _DEBLUR_THREAD_

#include "wx/progdlg.h"
#include "wx/thread.h"
#include "wx/event.h"
#include "wx/app.h"
#include "libdeblur.h"

using namespace cv;
using namespace std;
using namespace Eigen;

/*-------------------------------------------------------
 * class
 -------------------------------------------------------*/
class DeblurThread : public BlindDeblur, public wxThread {
  
 public:
 DeblurThread() : wxThread(wxTHREAD_DETACHED) {};
  virtual void *Entry();
  ~DeblurThread() {};
  
  void setCount(const int c);
  virtual bool checkInterrupt();
  
 private:
  int Count;  
};

#endif
