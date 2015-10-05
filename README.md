# Blind-Deblurring

This is a blind deblurring tool. It only works on Linux, and requires the following libraries:

 - [wxWidgets 3.0]
 - [OpenCV 2.4]
 - [fftw3.4]
 - [Eigen 3]

This tool is composed of a deblurring engine and GUI. Its engine is provided as a static library because source code of this engine is a top secret :-) Honestly speaking, it has huge numbers of bugs, so I have no plan to provide it.

I provide a binary image, so you can try it.

    $ ./blindDeblur

My goal is to make the best deblurring tool, and I thought it is easy  because it's just to make a solver of classical nonlinear optimization problem. But I had misunderstood. The destination is far away.


This tool is the first small step for me, and I'm going to improve it.


Screenshot
----------
**Original blurred image**
![alt text](http://www.interdb.jp/screenshot01.jpg)


**deblurring...**
![alt text](http://www.interdb.jp/screenshot2.jpg)

**after deblurring**
![alt text](http://www.interdb.jp/screenshot3.jpg)



[wxWidgets 3.0]: https://www.wxwidgets.org/
[OpenCV 2.4]: http://opencv.org/
[fftw3.4]: http://www.fftw.org/
[Eigen 3]: http://eigen.tuxfamily.org/index.php?title=Main_Page

