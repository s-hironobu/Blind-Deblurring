#include "deblur_thread.h"

void *DeblurThread::Entry() {
  int max = 100;
  try {
    exe(deblur->getImage());
  }
  catch (int err) {
    return NULL;
  }
  return NULL;
}

void DeblurThread::setCount(const int c) {
  /* dummy */ Count = c;
};

bool DeblurThread::checkInterrupt() { 
  return TestDestroy();
};
  


