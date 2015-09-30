#ifndef _CONFIG_
#define _CONFIG_

// default tool size
#define DEFAULT_WIDGET_WIDTH 800
#define DEFAULT_WIDGET_HEIGHT 700

// limitations of Image
#define MAX_FILE_SIZE 4*1024*1024 /* 4Mbyte */
#define MIN_IMAGE_H_SIZE 100
#define MIN_IMAGE_W_SIZE 100
#define MAX_IMAGE_H_SIZE 2000
#define MAX_IMAGE_W_SIZE 2000

// params for kernel display
#define LENGTH 13
#define MARGIN 20
#define SCALE 5

#define PSF_SIZE 32
#define DISP_PSF_SIZE 65

#define INIT_MODE 0

// default Image
#define DEFAULT_IMAGE_PATH "./sample-high.jpg"
#define DEFAULT_IMAGE "sample-high.jpg"

#endif
