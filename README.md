# Blind-Deblurring System


This is a blind deblurring system which is based on this paper: ["Blind deconvolution using alternating maximum a posteriori estimation with heavy-tailed priors"](https://users.soe.ucsc.edu/~milanfar/publications/conf/CAIP_paper_154.pdf).


I provide a [vagrant-box](https://atlas.hashicorp.com/s-hironobu/boxes/centos7-blind-deblur) that contains both binary and source code, so you can easily try it.

The most important remaining work is **to reduce ringing artifacts**.
Honestly speaking, I had been trying it until this summer when I had time, but I could not do it.
I welcome the challengers of this task :-)


In addition, the license of this system is *GPL2*.


## Description

This system runs on only Linux, and requires the following libraries:

 - [wxWidgets 3.0]
 - [OpenCV 2.4]
 - [fftw3.4]
 - [Eigen 3]

This system is composed of a deblurring engine and a GUI sub-system.

The engine is provided as a static library under the `lib` subdirectory, 
and the source code of it is also provided under the `libsrc` subdirectory.
If you recreate the static library, you have to execute `make` command on the *libsrc* subdirectory.


## Vagrant box

### Requirement

* [Vagrant](https://www.vagrantup.com/) 
* [VirtualBox](https://www.virtualbox.org/)
* vnc client

### How to use

#### [1] Create Vagrantfile

Create a Vargantfile as shown below:

```
# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
  config.vm.box = "s-hironobu/centos7-blind-deblur"
  config.vm.network "private_network", ip: "192.168.33.95"
  config.vm.network :forwarded_port, guest: 5901, host: 5901
end
```
#### [2] Start vm box

Run `vagrant up`.

```
$ vagrant up
```

#### [3] Access VM via vnc client

If you use OSX, run the following command:

```
$ open vnc://localhost:5901
```

Password is `vagrant`.

![Figure 1:](http://www.interdb.jp/blinddeblurring/blind-deblurring.png)

If you use other OS, run your vnc-client. 


#### [4] Open terminal and Run

Open a terminal and change directory to `blind-deblurring`.
Then, run `blindDeblur`.

```
[vagrant@localhost ~]$ cd blind-deblurring
[vagrant@localhost blind-deblurring]$  ./blindDeblur
```

Push `Deblur` button, and then adjust the sharpness of the deblurred image.


## Screenshot

**Original blurred image**
![alt text](http://www.interdb.jp/screenshot01.jpg)


**deblurring...**
![alt text](http://www.interdb.jp/screenshot2.jpg)

**Deblurred image**
![alt text](http://www.interdb.jp/screenshot3.jpg)


[wxWidgets 3.0]: https://www.wxwidgets.org/
[OpenCV 2.4]: http://opencv.org/
[fftw3.4]: http://www.fftw.org/
[Eigen 3]: http://eigen.tuxfamily.org/index.php?title=Main_Page

