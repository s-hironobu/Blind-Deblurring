CC = g++ -O2 -m64 -Wall
LIBS =  `pkg-config --libs opencv fftw3`
CFLAGS= `pkg-config --cflags opencv` -I /usr/local/include/eigen3/ -I /usr/local/include -I . -DEIGEN_NO_DEBUG

SLIB = libdeblur.a
OBJS = mat_utils.o edgetaper.o fft.o blind_deblur.o deblur.o
PROG = $(SRC:%.cpp=%) $(OBJS)
LIBDIR = ../lib

all: $(PROG)
	ar rvs $(LIBDIR)/$(SLIB) $(OBJS)
	strip -x $(LIBDIR)/$(SLIB)

##
## Blind Deblurring
##
blind_deblur.o: blind_deblur.cpp blind_deblur.h
	$(CC) $(CFLAGS) -c blind_deblur.cpp

deblur.o: deblur.cpp deblur.h
	$(CC) $(CFLAGS) -c deblur.cpp

mat_utils.o: mat_utils.c mat_utils.h
	$(CC) $(CFLAGS) -c mat_utils.c

edgetaper.o: edgetaper.c edgetaper.h
	$(CC) $(CFLAGS) -c edgetaper.c

fft.o: fft.c fft.h
	$(CC) $(CFLAGS) -c fft.c

clean: 
	rm -f *o *~ $(PROG) $(SLIB)
