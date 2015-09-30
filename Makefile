### Variables: ###
WX_CONFIG = wx-config
TOOLKIT = GTK

###
WX_RELEASE = `$(WX_CONFIG) --query-version`
WX_VERSION = `$(WX_CONFIG) --version`
$(WX_CONFIG) --toolkit=$(TOOLKIT)

wx_top_builddir = /usr/local/src/wxWidgets-3.0.2
#LIBDIRNAME = $(wx_top_builddir)/lib
LIBDIRNAME = /usr/local/lib

###
CPP = g++
CPPWARNINGS = -Wall -Wundef -Wunused-parameter -Wno-ctor-dtor-privacy -Woverloaded-virtual
LIBS = -lz -ldl -lm -llzma  `pkg-config --libs opencv gtk+-2.0 fftw3` -pthread ./lib/libdeblur.a
CPPFLAGS = -D_FILE_OFFSET_BITS=64 -DWX_PRECOMP -D__WX$(TOOLKIT)__ -O2 -fno-strict-aliasing \
	`pkg-config --cflags gtk+-2.0` `$(WX_CONFIG) --cppflags` \
	-I . -I/usr/local/include/eigen3

srcdir = ./src
objdir = ./obj
BLINDDEBLUR_OBJECTS = $(objdir)/deblur_thread.o $(objdir)/blindDeblur.o

PORTNAME = `$(WX_CONFIG) --query-toolkit`
WXUNICODEFLAG = u
COND_MONOLITHIC_0___WXLIB_RIBBON_p = \
	-lwx_$(PORTNAME)$(WXUNIVNAME)$(WXUNICODEFLAG)$(WXDEBUGFLAG)_ribbon-$(WX_RELEASE)
__WXLIB_RIBBON_p = $(COND_MONOLITHIC_0___WXLIB_RIBBON_p)

### Targets: ###
all: blindDeblur

blindDeblur: $(BLINDDEBLUR_OBJECTS) $(srcdir)/libdeblur.h $(srcdir)/config.h $(srcdir)/bitmaps.h
	$(CPP) -o $@ $(BLINDDEBLUR_OBJECTS) -L$(LIBDIRNAME) $(LIBS) $(__WXLIB_RIBBON_p) $(PLUGIN_ADV_EXTRALIBS)

$(objdir)/blindDeblur.o: $(srcdir)/blindDeblur.cpp
	$(CPP) -c -o $@ $(CPPFLAGS) $(srcdir)/blindDeblur.cpp

$(objdir)/deblur_thread.o: $(srcdir)/deblur_thread.cpp $(srcdir)/deblur_thread.h
	$(CPP) -c -o $@ $(CPPFLAGS) $(srcdir)/deblur_thread.cpp

clean: 
	rm -f ./*.o ./*~ $(srcdir)/*~ $(objdir)/*.o
	rm -f blindDeblur
