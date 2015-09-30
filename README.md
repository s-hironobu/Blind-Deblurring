# Blind-Deblurring

This is a blind deblurring tool. It only works on Linux, and requires the following libraries:

 - [wxWidgets 3.0]
 - [OpenCV 2.4]
 - [fftw3.4]
 - [Eigen 3]


I provide the source code of a part of GUI, and a static library of deblurring engine.

I have no plan to provide the source code of the deblurring engine since it is the top secret :-) Honestly speaking, the engine has huge numbers of bugs, so I can't release it yet.

My final goal is to make the best deblurring tool, and I thought it is easy  because it's just to make a solver of a classical nonlinear optimization program. But, as you can see, I had misunderstood.  The destination is far away.


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

