/*
 * Blind deblur tool
 *
 * Copyright (C) 2015  suzuki hironobu
 */

// ----------------------------------------------------------------------------
// declarations
// ----------------------------------------------------------------------------

/*
 * headers
 */
#include "wx/wx.h"
#include "wx/notifmsg.h"
#include "wx/progdlg.h"
#include "wx/file.h"
#include "wx/tokenzr.h"
#include "wx/combobox.h"

#include "config.h"
#include "libdeblur.h"
#include "deblur_thread.h"

/*
 * resources
 */

#include "bitmaps.h"

/*
 * namespace
 */
using namespace cv;
using namespace std;

/*
 *constants
 */
enum
{
  // menu items
  ID_menuFileLoadImage,
  ID_menuFileSaveImage,
  ID_menuFileQuit = wxID_EXIT,
  ID_menuAbout = wxID_ABOUT,

  // tool bar
  ID_toolBarLoadImage,
  ID_toolBarSaveImage,
  ID_toolBarHelp,
  ID_toolBarZoomIn,
  ID_toolBarZoomOut,
  ID_toolBarFitImage,
  ID_toolBarResizeImage,

  // operation panel
  ID_opFileLoadImage  = wxID_HIGHEST + 1,
  ID_opFileSaveImage,
  ID_opDeblurMode,
  ID_opSliderSharpness,
  ID_opCurTextSharpness,
#ifdef _BRIGHT_
  ID_opSliderBrightness,
  ID_opCurTextBrightness,
#endif
  ID_opDeblurImage,
  ID_opDeblurStop
};

enum 
  {
    DISP_ORIG,
    DISP_DEBLURED,
    DISP_FILTERED
  };

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

/*
 * store three kinds of images (original, deblured, filtered), and kernel data
 */
class ImageData : public wxObject
{
public:
  ImageData() {
    height = width = 0;
    isDemo = true;
    isChanged = false;
    deblur = NULL;  

    filepath.Clear();
    file.Clear();

    image_orig = new wxImage();
    image_deblured = new wxImage(); 
    image_filtered = new wxImage();
    image_kernel = new wxImage();
  };

  ~ImageData() {
    delete(deblur);
    filepath.Clear(); 
    file.Clear();
    delete(image_orig);
    delete(image_deblured);
    delete(image_filtered);
    delete(image_kernel);
  };

  bool SetImageData(wxString filepathname, wxString filename, wxBitmapType format, 
		    wxString &msg, bool isDemo = false);
  bool isSaved()  {
    return ((!isDemo && isChanged) ? true:false);
  }
  wxImage* GetImage() {
    return ((image_orig->IsOk())? image_orig : NULL);
  };
  wxImage* GetDebluredImage() {
    return ((image_deblured->IsOk())? image_deblured : NULL);
  };
  wxImage* GetFilteredImage() {
    return ((image_filtered->IsOk())? image_filtered : NULL);
  };
  Mat getKernelImage () {
    return deblur->getKernelImage();
  };

  bool GetImageSize(int *width, int *height);
  void DeblurData(const int mode, const float sharpness
#ifdef _BRIGHT_
		  , const float brightness
#endif
		  , const int psfSize
		  );  

  DeblurThread *GetDeblurThread() { return deblurThread;};

  void SharpImage (const float k);

  void setParam(const int mode) {
    deblur->setParam(mode);
  };

  wxString GetFileName() { return file;};
  wxString GetFilePath() { return filepath;};
  bool SaveImage(const wxString filepath) { return deblur->saveImage(filepath.mb_str());};

  // utils
  wxBitmapType GetFormat(wxString filename);
  bool CheckFileType(wxString filename);
  wxString GetSuffix(wxString filename);

private:
  Deblur *deblur;
  DeblurThread *deblurThread;

  wxString filepath;        // image's fileparh 
  wxString file;            // image's filename
  int height, width;        // image's height and width
  wxBitmapType imageFormat; // original image format, e.g, wxBITMAP_TYPE_JPEG, wxBITMAP_TYPE_GIF, etc.

  // 255x3channels (RGB)
  wxImage *image_orig;      // original image
  wxImage *image_deblured;  // deblured image
  wxImage *image_filtered;  // deblured+filtered image

  wxImage *image_kernel;    // kernel image

  bool isDemo, isChanged;
  
  void ClearData();
  void ReadDebluredImageFromMat(void);
  void ReadFilteredImageFromMat(void);
  void ReadKernelImageFromMat(void);
  void imageFromMat(Mat img, wxImage *image);
};

/*
 * Scrolled window
 */
class MyCanvas : public wxScrolledWindow
{
public:
  MyCanvas(wxFrame* parent, wxWindowID id, 
	   const wxPoint &pos, const wxSize &size);
  ~MyCanvas(){
    delete(imageData);
    dispImage.~wxBitmap();
  };

  void DisplayImage(void);
  void ResetScale(void) { Scale = 1.0;};

  //
  ImageData* GetImageData() const {return imageData;};

  // zoom operations
  void OnZoomIn(void);
  void OnZoomOut(void);
  void OnFitImage(void);
  void OnResizeImage(void);

private: 
  ImageData *imageData;

  float Scale;
  int width, height;
  wxBitmap dispImage;
  int dispMode;

  //   
  void paintEvent(wxPaintEvent & event);
  void render(wxDC& dc);

  // event handler
  void OnLDown(wxMouseEvent &WXUNUSED(event));
  void OnLUp(wxMouseEvent &WXUNUSED(event));

  wxDECLARE_EVENT_TABLE();
};

/*
 * main frame
 */
class MyFrame : public wxFrame
{
public:
  MyFrame(const wxString& title);

  ~MyFrame() {
    delete(m_canvas);
    delete(kernelImage);
    delete(kernel_window);
    delete(wxSliderSharpness);
#ifdef _BRIGHT_
    delete(wxSliderBrightness);
#endif
    delete(curTextSharpness);
#ifdef _BRIGHT_
    delete(curTextBrightness);
#endif
    currentDir.Clear();
    wildcards.Clear();
  };

  MyCanvas *GetCanvas() const {return m_canvas;};
  wxCriticalSection cs;
  int getPsfSize() { return psfSize;};

private:
  MyCanvas *m_canvas;

  wxImage *kernelImage;
  wxWindow *kernel_window;

  wxString progName;

  // param
  int psfSize;
  int mode;

  // radio box
  wxRadioBox *m_radioBox;

  // sliders
  wxSlider *wxSliderSharpness;
#ifdef _BRIGHT_
  wxSlider *wxSliderBrightness;
#endif

  int sliderValSharpness;
#ifdef _BRIGHT_
  int sliderValBrightness;
#endif

  wxTextCtrl *curTextSharpness;
#ifdef _BRIGHT_
  wxTextCtrl *curTextBrightness;
#endif
  // for fileDialog. 
  wxString currentDir;
  wxString wildcards;
  
  // create items
  void createMenu();
  void createToolBar();
  void initToolBar(wxToolBar* toolBar);
  void createControlPanel(wxBoxSizer *parent);
  void createCanvas(wxBoxSizer *parent);
  void createStatusBar();
  wxBoxSizer *createSlider(wxWindow *parent, wxSlider **slider, wxWindowID slider_id,
			   wxTextCtrl **textctrl, wxWindowID textctrl_id,
			   wxString label, int init, int min, int max, int *item);
  
  // convertors
  float convSharpness(const int sharpness);
#ifdef _BRIGHT_
  float convBrightness(const int brightness);
#endif

  // display kernel  
  void paintKernelEvent(wxPaintEvent &WXUNUSED(event));
  void displayKernelImage(void);
  void kernelRender(wxDC&  dc);
  void drawBitmap(wxDC&  dc);

  // event handlers
  void OnQuit(wxCommandEvent& event);
  void OnAbout(wxCommandEvent& event);
  void OnDeblurData(wxCommandEvent& event);
  void OnDeblurStop(wxCommandEvent& event);
  void OnLoadData(wxCommandEvent& event);
  void OnSaveData(wxCommandEvent& event);
  void OnZoomIn(wxCommandEvent &WXUNUSED(event));
  void OnZoomOut(wxCommandEvent &WXUNUSED(event));
  void OnFitImage(wxCommandEvent &WXUNUSED(event));
  void OnResizeImage(wxCommandEvent &WXUNUSED(event));

  void OnCheckRadioBox(wxCommandEvent& WXUNUSED(event));
  void OnUpdateSliderSharpness( wxCommandEvent &WXUNUSED(event));
#ifdef _BRIGHT_
  void OnUpdateSliderBrightness( wxCommandEvent &WXUNUSED(event));
#endif
  void OnUpdateCurTextSharpness(wxUpdateUIEvent& event);
#ifdef _BRIGHT_
  void OnUpdateCurTextBrightness(wxUpdateUIEvent& event);
#endif
  // others
  void postProduction();
  void setCurrentDir(wxString filepath, wxString filename);

  wxDECLARE_EVENT_TABLE();  
};

/*
 * top level application
 */
class MyApp : public wxApp
{
  MyFrame *main_frame;
public:
  virtual bool OnInit();
  MyFrame *GetFrame() const {return main_frame;};
};

// ----------------------------------------------------------------------------
// implement application
// ----------------------------------------------------------------------------

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
  if ( !wxApp::OnInit() )
    return false;
  wxInitAllImageHandlers();
  main_frame = new MyFrame("Blind deblurring tool");
  main_frame->Show(true);

  return true;
}

DECLARE_APP(MyApp)

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
/* menu
 */
 EVT_MENU(ID_menuFileLoadImage,  MyFrame::OnLoadData)
 EVT_MENU(ID_menuFileSaveImage,  MyFrame::OnSaveData)
 EVT_MENU(ID_menuFileQuit,  MyFrame::OnQuit)
 EVT_MENU(ID_menuAbout, MyFrame::OnAbout)
/* toolbar
 */
 EVT_MENU(ID_toolBarLoadImage, MyFrame::OnLoadData)
 EVT_MENU(ID_toolBarSaveImage, MyFrame::OnSaveData)
 EVT_MENU(ID_toolBarZoomIn, MyFrame::OnZoomIn)
 EVT_MENU(ID_toolBarZoomOut, MyFrame::OnZoomOut)
 EVT_MENU(ID_toolBarFitImage, MyFrame::OnFitImage)
 EVT_MENU(ID_toolBarResizeImage, MyFrame::OnResizeImage)
/* control panel
 */
// EVT_BUTTON(ID_opFileLoadImage, MyFrame::OnLoadData)
// EVT_BUTTON(ID_opFileSaveImage, MyFrame::OnSaveData)
 EVT_CHECKBOX(ID_opDeblurMode, MyFrame::OnCheckRadioBox)
 EVT_RADIOBOX(ID_opDeblurMode, MyFrame::OnCheckRadioBox)

 EVT_SLIDER(ID_opSliderSharpness, MyFrame::OnUpdateSliderSharpness)
#ifdef _BRIGHT_
 EVT_SLIDER(ID_opSliderBrightness, MyFrame::OnUpdateSliderBrightness)
#endif

 EVT_UPDATE_UI(ID_opCurTextSharpness, MyFrame::OnUpdateCurTextSharpness)
#ifdef _BRIGHT_
 EVT_UPDATE_UI(ID_opCurTextBrightness, MyFrame::OnUpdateCurTextBrightness)
#endif
 EVT_BUTTON(ID_opDeblurImage, MyFrame::OnDeblurData)
 EVT_BUTTON(ID_opDeblurStop, MyFrame::OnDeblurStop)
/* paint
 */
 EVT_PAINT(MyFrame::paintKernelEvent)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyCanvas, wxScrolledWindow)
 EVT_PAINT(MyCanvas::paintEvent)
 EVT_LEFT_DOWN(MyCanvas::OnLDown)
 EVT_LEFT_UP(MyCanvas::OnLUp)
wxEND_EVENT_TABLE()

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// class ImageData
// ----------------------------------------------------------------------------

/*
 * clear data
 */
void ImageData::ClearData()
{
  if(image_orig->IsOk())
    image_orig->Destroy();
  if(image_deblured->IsOk())
    image_deblured->Destroy();
  if(image_filtered->IsOk())
    image_filtered->Destroy();
  if(image_kernel->IsOk())
    image_kernel->Destroy();

  height = 0;
  width = 0;
  isDemo = false; // default
  isChanged = false;

  filepath.Clear();
  file.Clear();
}


/*
 *
 */
bool ImageData::GetImageSize(int *width, int *height)
{
  if (image_orig->IsOk()) {
    *width = image_orig->GetWidth();
    *height = image_orig->GetHeight();
    return true;
  }
  else{
    *width = -1;
    *height = -1;
  }
  return false;
}

/*
 *
 */
bool ImageData::SetImageData(wxString filepathname, wxString filename, 
			     wxBitmapType format, 
			     wxString &msg, bool Demo)
{
  // check file size
  wxFile *tmp = new wxFile(filepathname, wxFile::read);

  if (!tmp->Exists(filepathname)) {
    msg = wxString::Format("Error: %s not found.", filename.mb_str());
    return false;
  }
  wxFileOffset len = tmp->Length();
  delete(tmp);

  if (MAX_FILE_SIZE < len) {
    msg = wxString::Format("Error: This file is over %d byte", MAX_FILE_SIZE);
    return false;
  }

  // set new image
  ClearData(); // clear old image

  // TODO : only check size?
  image_orig->LoadFile(filepathname, format);

  width = image_orig->GetWidth();
  height = image_orig->GetHeight();

  if (width < MIN_IMAGE_W_SIZE || height < MIN_IMAGE_H_SIZE) {
    msg = wxString::Format("Error: This image is less than  (%d x %d)",
			   MIN_IMAGE_W_SIZE, MIN_IMAGE_H_SIZE);
    return false;
  }
  else if (MAX_IMAGE_W_SIZE < width || MAX_IMAGE_H_SIZE < height) {
    msg = wxString::Format("Error: This image is greater than  (%d x %d)",
			   MAX_IMAGE_W_SIZE, MAX_IMAGE_H_SIZE);
    return false;
  }
  
  // create deblur instance
  if (deblur != NULL)
    delete(deblur);
  
  try {
    deblur = new Deblur(filepathname.mb_str());
    deblur->initKernelImage(PSF_SIZE);
  }
  catch (char *err) {
    msg = wxString::Format("%s", err);
    return false;
  }

  //   
  ReadDebluredImageFromMat();
  ReadFilteredImageFromMat();
  
  // initialize params
  isDemo = Demo;
  isChanged = false;
  filepath = filepathname;
  file = filename;
  imageFormat = format;

  return true;
}

/*
 * convert from OpenCV's Mat(cv_8U) to wxImage
 */
void ImageData::imageFromMat(Mat img, wxImage *image)
{
  if (image->IsOk())
    image->Destroy();
  image->Create(width, height);

  unsigned char * pixels = image->GetData();

  // Assume that the Mat::img is 2dim, 3 channels, and CV_8U (8bit=255).
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      int index = (3 * x) + (3 * y * width);
      int key = (y * img.step) + (x * img.elemSize());

      pixels[index+0] = img.data[key + 2];
      pixels[index+1] = img.data[key + 1];
      pixels[index+2] = img.data[key + 0];
    }
  }

  isChanged = true;
}

/*
 *
 */
void ImageData::ReadDebluredImageFromMat(void)
{
  imageFromMat(deblur->getDebluredImage(), image_deblured);
}

/*
 *
 */
void ImageData::ReadFilteredImageFromMat(void)
{
  imageFromMat(deblur->getFilteredImage(), image_filtered);
}

/*
 *
 */
void ImageData::ReadKernelImageFromMat(void)
{
  imageFromMat(deblur->getKernelImage(), image_kernel);
}

/*
 *
 */
void ImageData::DeblurData( const int mode, const float sharpness
#ifdef _BRIGHT_
			   , const float brightness
#endif
			   , const int psfSize			   
			   )
{
  deblurThread = new DeblurThread();
  if (deblurThread->Create() != wxTHREAD_NO_ERROR) {
    wxLogError(wxT("Can't create thread!"));
  }
  // TODO  set virtual function::: 
  deblur->initKernelImage(psfSize);

  deblurThread->setDeblur(deblur);
  deblurThread->setMode(mode);
  deblurThread->setPsfSize(psfSize);
  deblurThread->setCount(0);
  deblurThread->Run();

  static const int max = 100;
  wxProgressDialog dialog(wxT("deblur dialog window"),
			  wxT(""),
                          max,    // range
                          wxGetApp().GetFrame(),   // parent
                          wxPD_CAN_ABORT |
                          //wxPD_CAN_SKIP |
                          wxPD_AUTO_HIDE |
                          wxPD_APP_MODAL |
                          wxPD_ELAPSED_TIME |
                          wxPD_ESTIMATED_TIME |
                          //wxPD_REMAINING_TIME |
                          wxPD_SMOOTH
                          );
  wxString msg;
  bool skip = false;
  msg = "Deblurring....\n";

  int count;
  {
    wxCriticalSectionLocker enter(wxGetApp().GetFrame()->cs);
    count = deblurThread->getCount();
  }

  while (count < 100) {
    {
      wxCriticalSectionLocker enter(wxGetApp().GetFrame()->cs);
      count = deblurThread->getCount();
    }

    bool flg =  dialog.Update(count, msg, &skip);
    wxMilliSleep(200);
    
    if (flg == false) {
      deblurThread->Delete();
      //      delete(deblurThread);
      break;
    }
  }

  deblur->SharpFilter(sharpness);
  
  ReadDebluredImageFromMat();
  ReadFilteredImageFromMat();
  ReadKernelImageFromMat();
}

/*
 *
 */
void ImageData::SharpImage(const float sharpness)
{
  deblur->SharpFilter(sharpness);
  ReadFilteredImageFromMat();
}


/*
 * util
 */
wxBitmapType ImageData::GetFormat(wxString filename)
{
  wxBitmapType ret = (wxBitmapType)NULL;
  wxStringTokenizer tkz(filename, wxT("."));
  wxString suffix;
  while (tkz.HasMoreTokens()) {
    suffix = tkz.GetNextToken();
  }

  if (wxStrcmp(suffix.Upper(), "BMP") == 0)
    ret = wxBITMAP_TYPE_BMP;
  else if (wxStrcmp(suffix.Upper(), "GIF") == 0)
    ret = wxBITMAP_TYPE_GIF;
  else if (wxStrcmp(suffix.Upper(), "JPG") == 0
	   || wxStrcmp(suffix.Upper(), "JPEG") == 0)
     ret = wxBITMAP_TYPE_JPEG;
  else if (wxStrcmp(suffix.Upper(), "PNG") == 0)
    ret = wxBITMAP_TYPE_PNG;
  else if (wxStrcmp(suffix.Upper(), "PCX") == 0)
    ret = wxBITMAP_TYPE_PCX;
  else if (wxStrcmp(suffix.Upper(), "PNM") == 0)
    ret = wxBITMAP_TYPE_PNM;
  else if (wxStrcmp(suffix.Upper(), "TIFF") == 0)
    ret= wxBITMAP_TYPE_TIFF;

  tkz.~wxObject();
  suffix.Clear();

  return ret;
}

bool ImageData::CheckFileType(wxString filename)
{
  if (GetFormat(filename) == (wxBitmapType)NULL)
    return false;
  return true;
}

wxString ImageData::GetSuffix(wxString filename)
{
  wxStringTokenizer tkz(filename, wxT("."));
  wxString suffix;
  while (tkz.HasMoreTokens()) {
    suffix = tkz.GetNextToken();
  }

  tkz.~wxObject();

  return suffix;
}

// ----------------------------------------------------------------------------
// class MyCanvas
// ----------------------------------------------------------------------------

/*
 * constructor
 */
MyCanvas::MyCanvas(wxFrame* parent, wxWindowID id, 
		   const wxPoint &pos, const wxSize &size) :
  wxScrolledWindow(parent, id, pos, size, wxSB_VERTICAL|wxSB_HORIZONTAL)
{
  wxString msg;
  imageData = new ImageData();

  if (imageData->SetImageData(wxT(DEFAULT_IMAGE_PATH), wxT(DEFAULT_IMAGE),
			      wxBITMAP_TYPE_JPEG, msg, true)) {
    imageData->GetImageSize(&width, &height);
    ResetScale();  //Scale = 1.0;
    dispMode = DISP_FILTERED;  
    SetScrollbars(1, 1, width*3, height*3, true, true);

    msg.Clear();
  }
  else {
    fprintf(stderr, "Fatal Error: Default Image data(%s) not found!\n", DEFAULT_IMAGE_PATH);
    exit(-1);
  }
}

/*
 * display initial image. It is only invoked at the first time.
 */ 
void MyCanvas::paintEvent(wxPaintEvent &WXUNUSED(event))
{
  wxClientDC dc(this);
  PrepareDC(dc);
  // Don't execute Refresh() here!
  render(dc);  
}

/*
 * render image
 */
void MyCanvas::render(wxDC&  dc)
{
  wxImage* image;
  switch (dispMode) 
    {
    case DISP_ORIG:
      image = imageData->GetImage();
      break;
    case DISP_DEBLURED:
      image = imageData->GetDebluredImage();
      break;
    default: // DISP_FILTERED
      image = imageData->GetFilteredImage();
    }
  
  dispImage = wxBitmap(image->Scale( width*Scale, height*Scale));
  dc.DrawBitmap( dispImage, 0, 0, false );
}

/*
 *
 */
void MyCanvas::DisplayImage(void)
{
  wxClientDC dc(this);
  PrepareDC(dc);
  Refresh();
  render(dc);  
}

/*
 * zoom in
 */
void MyCanvas::OnZoomIn(void)
{
  if (Scale < 2.5)
    Scale += 0.1;
  else
    Scale = 2.5;

  DisplayImage();
}

/*
 * zoom out
 */
void MyCanvas::OnZoomOut()
{
  if (0.3 < Scale) 
    Scale -= 0.1;
  else 
    Scale = 0.3;
  DisplayImage();
}

/*
 * fit display size
 */
void MyCanvas::OnFitImage()
{
  int imageWidth, imageHeight;
  int visibleWidth, visibleHeight;
  
  Scale = 1.0;
  GetSize(&visibleWidth, &visibleHeight);

  wxImage* image = imageData->GetImage();
  if (image->IsOk()) {
    float wratio = ((float)visibleWidth / (float)width);
    float hratio = ((float)visibleHeight / (float)height);
    Scale = (wratio < hratio) ? wratio : hratio;
    DisplayImage();
  }
}

/*
 * display original sized image
 */
void MyCanvas::OnResizeImage()
{
  Scale = 1.0;
  DisplayImage();
}

/*
 * display original image
 */
void MyCanvas::OnLDown(wxMouseEvent &WXUNUSED(event))
{
  dispMode = DISP_ORIG;
  DisplayImage();
  wxGetApp().GetFrame()->SetStatusText("Your original image.");
}

/*
 * display deblured+filtered image
 */
void MyCanvas::OnLUp(wxMouseEvent &WXUNUSED(event))
{
  dispMode = DISP_FILTERED;
  DisplayImage();
  wxGetApp().GetFrame()->SetStatusText("Deblured image.");
}


// ----------------------------------------------------------------------------
// class MyFrame
// ----------------------------------------------------------------------------

/*
 * constructor
 */
MyFrame::MyFrame(const wxString& title)
  : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition,
	    wxSize(DEFAULT_WIDGET_WIDTH, DEFAULT_WIDGET_HEIGHT), wxDEFAULT_FRAME_STYLE)
{
  // set the program name
  progName = wxT("MagicDeblur");

  // init params
  psfSize = PSF_SIZE;
  mode = INIT_MODE;

  // set the frame icon
  SetIcon(wxICON(sample));
  
  // create menu
  createMenu();

  // create toolbar
  createToolBar();

  // main area
  wxBoxSizer *topSizer = new wxBoxSizer(wxVERTICAL);
  createControlPanel(topSizer);
  createCanvas(topSizer);
  SetSizer(topSizer);
  
  // status bar
  createStatusBar();

  // for display kernel image
  kernelImage = new wxImage();
  displayKernelImage();

  // default directory for fileDialog
  currentDir = wxT("./");
  wildcards = wxString::Format(wxT("JPEG file (*jpeg;*jpg)|*.jpg;*jpeg|\
BMP or GIF file (*.bmp;*.gif)|*.bmp;*.gif|			       \
PNG file (*.png)|*.png"),
			       wxFileSelectorDefaultWildcardStr,
			       wxFileSelectorDefaultWildcardStr,
			       wxFileSelectorDefaultWildcardStr
			       );  
}

/*-------------------------------
 * create items
 --------------------------------*/
/*
 * create menu
 */
void MyFrame::createMenu()
{
  // create File
  wxMenu *menuFile = new wxMenu;
  menuFile->Append(ID_menuFileLoadImage, wxT("&Load image"));
  menuFile->Append(ID_menuFileSaveImage, wxT("&Save image"));
  menuFile->AppendSeparator();
  menuFile->Append(ID_menuFileQuit, wxT("&Quit"));
  
  // create Help
  wxMenu *menuHelp = new wxMenu;
  menuHelp->Append(ID_menuAbout, "&About", "Show about dialog");
  
  // create MenuBar
  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append(menuFile, "&File");
  menuBar->Append(menuHelp, "&Help");

  SetMenuBar(menuBar);
}

/*
 * create toolbar
 */
void MyFrame::createToolBar()
{
  CreateToolBar(wxNO_BORDER | wxTB_FLAT | wxTB_HORIZONTAL);
  initToolBar(GetToolBar());
}

/*
 * initialize tool bar
 */
void MyFrame::initToolBar(wxToolBar* toolBar)
{
    const int maxBitmaps = 9;
    wxBitmap* bitmaps[maxBitmaps];

    bitmaps[0] = new wxBitmap( fileopen_xpm );
    bitmaps[1] = new wxBitmap( filesave_xpm );
    bitmaps[2] = new wxBitmap( help_xpm );
    bitmaps[3] = new wxBitmap( position_left_xpm );
    bitmaps[4] = new wxBitmap( position_top_xpm );
    bitmaps[5] = new wxBitmap( zoom_in_xpm );
    bitmaps[6] = new wxBitmap( zoom_out_xpm );
    bitmaps[7] = new wxBitmap( zoom_fit_best_xpm );
    bitmaps[8] = new wxBitmap( zoom_original_xpm );

    toolBar->AddTool(ID_toolBarLoadImage, wxEmptyString, *(bitmaps[0]), wxS("Load new file"));
    toolBar->AddTool(ID_toolBarSaveImage, wxEmptyString, *bitmaps[1], wxS("Save image"));
    toolBar->AddSeparator();
    toolBar->AddTool(ID_toolBarZoomIn, wxEmptyString, *bitmaps[5], wxS("Zoom in"));
    toolBar->AddTool(ID_toolBarZoomOut, wxEmptyString, *bitmaps[6], wxS("Zoom out"));
    toolBar->AddSeparator();
    toolBar->AddTool(ID_toolBarFitImage, wxEmptyString, *bitmaps[7], wxS("Fitting"));
    toolBar->AddTool(ID_toolBarResizeImage, wxEmptyString, *bitmaps[8], wxS("Original Size"));

    toolBar->Realize();
    int i;
    for (i = 0; i < maxBitmaps; i++)
      delete(bitmaps[i]);
}

/*
 * create control panel
 */
void MyFrame::createControlPanel(wxBoxSizer *parent)
{
  wxBoxSizer *topSizer = new wxBoxSizer(wxHORIZONTAL);

#ifdef _TEST_
  // load/save file 
  wxStaticBox *file_box = new wxStaticBox(this, wxID_ANY, wxT("File"));  
  wxStaticBoxSizer *file_sizer = new wxStaticBoxSizer(file_box, wxVERTICAL);
  {
    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
    wxButton *button1 = new wxButton(file_box, ID_opFileLoadImage, wxT("Load"));
    wxButton *button2 = new wxButton(file_box, ID_opFileSaveImage, wxT("Save"));
    sizer->Add(button1);
    sizer->Add(button2);
    file_sizer->Add(sizer);
  }
  topSizer->Add(file_sizer, 0, wxALL, 15);
#endif

  // Blind Deconvolution
  wxStaticBox *psf_box = new wxStaticBox(this, wxID_ANY, wxT("Deblurring"));
  wxStaticBoxSizer *psf_sizer = new wxStaticBoxSizer(psf_box, wxHORIZONTAL);
  
  {
    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);

    static const wxString deblurMode[] = {
      wxT("High"),
      wxT("Medium"),
      wxT("Low")
    };
    m_radioBox = new wxRadioBox(psf_box, ID_opDeblurMode, wxT("Deblurring power:"),
				wxDefaultPosition, wxDefaultSize,
				WXSIZEOF(deblurMode), deblurMode,
				1, wxRA_SPECIFY_COLS);
   sizer->Add(m_radioBox);
   psf_sizer->Add(sizer);
  }


  {
    // display kernel
    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
    kernel_window = new wxWindow(psf_box, wxID_ANY, wxDefaultPosition,
				 wxSize(DISP_PSF_SIZE+MARGIN,DISP_PSF_SIZE+MARGIN));
    sizer->Add(kernel_window, 0, wxALIGN_RIGHT|wxALL, 5);
    psf_sizer->Add(sizer);
  }  

  {
    // exec button
    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
    wxButton *button = new wxButton(psf_box, ID_opDeblurImage, wxT("Deblur"));
    sizer->Add(button);
    psf_sizer->Add(sizer);
  }
  topSizer->Add(psf_sizer, 0, wxALL, 15);

  // filter
  wxStaticBox *filter_box = new wxStaticBox(this, wxID_ANY, wxT("Post-Production"));
  wxStaticBoxSizer *filter_sizer = new wxStaticBoxSizer(filter_box, wxVERTICAL);
  {
    wxBoxSizer *sl_sizer1 = createSlider(this, &wxSliderSharpness, ID_opSliderSharpness, 
					 &curTextSharpness,
					 ID_opCurTextSharpness, wxT("sharpness"),
					 0, 0, 30, &sliderValSharpness);
    filter_sizer->Add(sl_sizer1);

#ifdef _BRIGHT_
    wxBoxSizer *sl_sizer2 = createSlider(this, &wxSliderBrightness, ID_opSliderBrightness,
					 &curTextBrightness,
					 ID_opCurTextBrightness, wxT("brightness"),
					 20, 0, 100, &sliderValBrightness);
    filter_sizer->Add(sl_sizer2);
#endif
  }
  topSizer->Add(filter_sizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 15);

  // add topSizer
  parent->Add(topSizer);
}

/*
 * create image display canvas
 */
void MyFrame::createCanvas(wxBoxSizer *topSizer)
{
  m_canvas = new MyCanvas(this, wxID_ANY, wxDefaultPosition, wxSize(10,10));
  topSizer->Add(m_canvas, 1, wxALL | wxEXPAND);
}

/*
 * create status bar
 */
void MyFrame::createStatusBar()
{
  CreateStatusBar(2);
  SetStatusText("This is a Blind deblurring tool.");
}

/*
 * create slider
 */
wxBoxSizer*  MyFrame::createSlider(wxWindow *parent, wxSlider **slider, wxWindowID slider_id,
                                   wxTextCtrl **textctrl, wxWindowID textctrl_id,
                                   wxString label, int init, int min, int max, int *item)
{
  wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  
  wxStaticText *sg_text = new wxStaticText(parent, wxID_STATIC, label, wxDefaultPosition, 
					   wxSize(80, -1), 0);
  sizer->Add(sg_text, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxALL, 5);

  *textctrl = new wxTextCtrl(parent, (wxWindowID)textctrl_id, wxEmptyString, 
			     wxDefaultPosition, wxSize(10,-1), wxTE_READONLY);
  sizer->Add(*textctrl, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxALL,5);

  *slider = new wxSlider(parent, slider_id, init, min, max, wxDefaultPosition, 
			 wxSize(100, wxDefaultCoord), wxSL_AUTOTICKS);
  *item = init; // set initial value

  sizer->Add(*slider, 0, wxALIGN_LEFT|wxALL, 5);
  
  return sizer;
}

/*-------------------------------
 * convertors
 --------------------------------*/
float MyFrame::convSharpness(const int sharpness)
{
  return (float)(sharpness / (float)10.0);
}

#ifdef _BRIGHT_
float MyFrame::convBrightness(const int brightness)
{
  return (float)(brightness / (float)10.0);
}
#endif


/*-------------------------------
 * display kernel
 --------------------------------*/
void MyFrame::paintKernelEvent(wxPaintEvent &WXUNUSED(event))
{
  wxClientDC dc(kernel_window);
  PrepareDC(dc);
  drawBitmap(dc);
}

void MyFrame::displayKernelImage(void)
{
  wxClientDC dc(kernel_window);
  PrepareDC(dc);
  Refresh();
  kernelRender(dc);
}

void MyFrame::kernelRender(wxDC&  dc)
{
  MyCanvas *m_canvas = NULL;
  m_canvas = this->GetCanvas();
  ImageData *imageData = NULL;
  imageData = m_canvas->GetImageData();
  Mat kernel = imageData->getKernelImage();

  if(kernelImage->IsOk())
    kernelImage->Destroy();

  kernelImage->Create(psfSize, psfSize);

  assert(kernel.rows == psfSize && kernel.cols == psfSize);

    unsigned char *pixels = kernelImage->GetData();

  for (int x = 0; x < psfSize; x++) {
    for (int y = 0; y < psfSize; y++) {
      int key = 3*x + 3*psfSize*y;
      unsigned char val =  kernel.data[y * kernel.step + x * kernel.elemSize() + 0];
      pixels[key+0] = val;
      pixels[key+1] = val;
      pixels[key+2] = val;
    }
  }

  drawBitmap(dc);
}

/*
 *
 */
void MyFrame::drawBitmap(wxDC&  dc)
{
  dc.DrawBitmap(wxBitmap(kernelImage->Rescale( DISP_PSF_SIZE, DISP_PSF_SIZE, wxIMAGE_QUALITY_NEAREST)),
		MARGIN/2, MARGIN/2, false );
}

/*-------------------------------
 * event handlers
 --------------------------------*/
/*
 * quit this program
 */
void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
  Close(true);
}

/*
 * create about dialog
 */
void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
  wxMessageBox(wxString::Format
	       (
		"Welcome to %s!\n"
		"\n"
		"This is a Blind-Deblurring Tool\n"
		"running under %s.",
		progName,
		wxGetOsDescription()
		),
	       "About",
	       wxOK | wxICON_INFORMATION,
	       this);
}

/*
 *
 */
void MyFrame::OnSaveData(wxCommandEvent& event)
{
  ImageData *imageData = this->GetCanvas()->GetImageData();

  wxString errmsg;
  wxString filename = wxString::Format(wxT("Untitled.%s"), 
				       imageData->GetSuffix(imageData->GetFileName()));
  // TODO message:
  wxFileDialog dialog(this, wxT("Save the image"),
		      currentDir, filename, wildcards,		      
		      wxFD_SAVE|wxFD_OVERWRITE_PROMPT);

  try {    
    if (dialog.ShowModal() == wxID_OK)
      {
	if (imageData->CheckFileType(dialog.GetPath())) {
	  if(imageData->SaveImage(dialog.GetPath())) {
	    setCurrentDir(dialog.GetPath(), dialog.GetFilename());
	  }
	  else {
	    errmsg = "Error: Couldn't save the image.";
	    throw errmsg;
	  }	  
	}
	else {
	  errmsg = wxString::Format("Error: \'%s\' is an invalid filename.", dialog.GetFilename());
	  throw errmsg;
	}
      }
    else {
      // Cancel button is pressed.
      SetStatusText("Cancelled File save operation.");     
    }
    errmsg.Clear();
    filename.Clear();
    
    event.Skip(true);
  }
  catch (wxString err) {
    wxLogMessage(wxT("%s"), err.mb_str());
    SetStatusText("Error: Couldn't save the image.");
    if (event.GetId() == ID_toolBarLoadImage
	|| event.GetId() == ID_menuFileLoadImage)
      throw err; 
  }
}

/*
 *
 */
void MyFrame::OnLoadData(wxCommandEvent& event)
{
  MyCanvas *m_canvas = this->GetCanvas();
  ImageData *imageData = m_canvas->GetImageData();

  try {
    // Is current image data saved?
    if(imageData->isSaved()) {
      wxMessageDialog dialog(this, "Save the current image before opening a new one?\n",
			     "Open new image", wxCENTER | wxNO_DEFAULT | wxYES_NO );

      dialog.SetYesNoLabels("&Yes", "Don't Save");
      
      wxString extmsg = "If you don't save the image, it will be lost.\n";
      dialog.SetExtendedMessage(extmsg);
      extmsg.Clear();
      
      switch ( dialog.ShowModal() )  {
      case wxID_YES:
	OnSaveData(event);
	break;
      case wxID_NO:
      default:
	;      
      }    
    }
    
    // create file open dialog
    wxFileDialog dialog(this, wxT("Testing open a file dialog"),
			currentDir, wxEmptyString, wildcards, wxFD_OPEN);
    // load image and display it.
    if (dialog.ShowModal() == wxID_OK) {

      wxString msg;
      if (imageData->SetImageData(dialog.GetPath(), dialog.GetFilename(),
				  imageData->GetFormat(dialog.GetPath()), msg, false)) {
	wxLogStatus(wxT("Load %s\n"), dialog.GetPath());
	m_canvas->ResetScale();      
	m_canvas->DisplayImage();
	
	setCurrentDir(dialog.GetPath(), dialog.GetFilename());
      }
      else {
	wxNotificationMessage n("Error", wxString::Format("%s", msg));
	n.SetFlags(wxICON_ERROR);
	n.Show(4); // Just for testing, use 4 second delay.      
      }
      msg.Clear();
    }
    event.Skip(false);

    SetStatusText("Image is saved");
    
  } catch(wxString err) {
    // terminate this method.
    SetStatusText("Error: Image could not saved.");
  }
}

/*
 *
 */
void MyFrame::OnDeblurStop(wxCommandEvent& WXUNUSED(event))
{
  ImageData *m_image = this->GetCanvas()->GetImageData();
  DeblurThread *deblurThread = m_image->GetDeblurThread();
  deblurThread->Delete();
  delete(deblurThread);
}

/*
 *
 */
void MyFrame::OnCheckRadioBox(wxCommandEvent& WXUNUSED(event)) {
    mode = m_radioBox->GetSelection();
}

/*
 *
 */
void MyFrame::OnDeblurData(wxCommandEvent& WXUNUSED(event))
{
  SetStatusText("Deblurring....");

  MyCanvas *m_canvas = this->GetCanvas();
  m_canvas->GetImageData()->DeblurData(mode,
				       convSharpness(sliderValSharpness)
#ifdef _BRIGHT_
				       ,convBrightness(sliderValBrightness)
#endif
				       ,psfSize
				       );
  // Display Deblured Image on the myCanvas
  m_canvas->DisplayImage();
  // Display Kernel on the myFrame
  displayKernelImage();

  SetStatusText("Deblurring....Finished.");
}

/*
 *
 */
void MyFrame::OnZoomIn(wxCommandEvent &WXUNUSED(event))
{
  this->GetCanvas()->OnZoomIn();
}

void MyFrame::OnZoomOut(wxCommandEvent &WXUNUSED(event))
{
  this->GetCanvas()->OnZoomOut();
}

void MyFrame::OnFitImage(wxCommandEvent &WXUNUSED(event))
{
  this->GetCanvas()->OnFitImage();
}

void MyFrame::OnResizeImage(wxCommandEvent &WXUNUSED(event))
{
  this->GetCanvas()->OnResizeImage();
}

void MyFrame::OnUpdateSliderSharpness( wxCommandEvent &WXUNUSED(event) )
{
  sliderValSharpness = wxSliderSharpness->GetValue();
  postProduction();
}

#ifdef _BRIGHT_
void MyFrame::OnUpdateSliderBrightness( wxCommandEvent &WXUNUSED(event) )
{
  sliderValBrightness = wxSliderBrightness->GetValue();
  postProduction();
}
#endif

void MyFrame::OnUpdateCurTextSharpness(wxUpdateUIEvent& event)
{    
  event.SetText( wxString::Format(wxT("%.1f"),
				  convSharpness(wxSliderSharpness->GetValue())));
  event.Skip(false);
}

#ifdef _BRIGHT_
void MyFrame::OnUpdateCurTextBrightness(wxUpdateUIEvent& event)
{    
  event.SetText( wxString::Format(wxT("%.1f"),
				  convBrightness(wxSliderBrightness->GetValue())));
  event.Skip(false);
}
#endif


/*-------------------------------
 * others
 --------------------------------*/
void MyFrame::postProduction()
{
  MyCanvas *m_canvas = this->GetCanvas();
  ImageData *imageData = m_canvas->GetImageData();
  
  imageData->SharpImage(convSharpness(sliderValSharpness));
  m_canvas->DisplayImage();
}

void MyFrame::setCurrentDir(wxString filepath, wxString filename)
{
  currentDir = filepath;
  currentDir.RemoveLast(filename.Len());
}    

// EOF
