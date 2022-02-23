# High-quality Ellipse Detection
## 1. Illustration
- This is the source code for the paper [Arc-support Line Segments Revisited: An Efficient and High-quality Ellipse Detection](https://arxiv.org/abs/1810.03243). ***Important: Please use the citation of our IEEE TIP version instead of arXiv version***.
- The main contribution of the proposed ellipse detector is to both accurately and efficiently detect ellipses in images, which is universally considered as a tough and long-standing problem in ellipse detection field. The proposed ellipse detector owns the features of *high localization accuracy, efficiency, robustness*, and *stability*, which comprehensively yields high-quality ellipse detection performance in front of real-world images. 
- There are only *two* extrinsic parameters, namely the elliptic angular coverage $T_{ac}$ and the ratio of support inliers $T_{r}$, which enables the proposed ellipse detector to be conveniently used and applied in real applications. In addition, the *specified_polarity* option can help users find the polarity-specific ellipses in the image. The default parameters $T_{ac} = 165^o$ and $T_{r} = 0.6$ are used for comparison experiments in our paper.  
- The source code is free for academic use. Please cite our paper if you use the source code, thanks.

## 2. Requirements
- MATLAB
- OpenCV (Version 2.4.9)
- 64-bit Windows Operating System


## 3. How to use
- Firstly, compile the file "generateEllipseCandidates.cpp" in MATLAB on your computer to generate the mex file "generateEllipseCandidates.mexw64" with the following command:  
  
  ---
  mex generateEllipseCandidates.cpp -IF:\OpenCV\opencv2.4.9\build\include -IF:\OpenCV\opencv2.4.9\build\include\opencv -IF:\OpenCV\opencv2.4.9\build\include\opencv2 -LF:\OpenCV\opencv2.4.9\build\x64\vc11\lib -IF:\Matlab\settlein\extern\include -LF:\Matlab\settlein\extern\lib\win64\microsoft -lopencv_core249 -lopencv_highgui249 -lopencv_imgproc249 -llibmwlapack.lib  
  
  ---
  Notably, the corresponding software paths of OpenCV and MATLAB, namely the "F:\OpenCV\opencv2.4.9\" and "F:\Matlab\settlein\", should be replaced to your own.  
- Secondly, run the demo file "LCS_ellipse.m".

## 3.5 How to use (extended)
- Before ellipse detection can be performed, the "generateEllipseCandidates.cpp" file needs to be compiled in order to generate a mex file. On Windows, this will create a .mexw64 file, and on macOS, it will create a .mexmaci64 file. Linux functionaility has not been tested.

  First, mex needs to know what compiler to use. On macOS, this will likely be with clang through XCode, and on Windows can be done with either mingW or one a compiler from Visual Studio. There have been reported issues with attempting to use mingW and some versions of Visual Studio, but the 2012 version appears to have the most success.

---
    mex -setup C++
    
---

- The command to mex the C++ code is below. While it is recommended to do this in the MATLAB terminal, it can be done in any terminal so long as MATLAB is correctly added to the PATH:

---
    mex -v generateEllipseCandidates.cpp -I'/path_to_opencv2_version/include' -I'/path_to_opencv2_version/include/opencv' -I'/path_to_opencv2_version/include/opencv2' -L'/path_to_opencv2_version/lib' -I'/path_to_MATLAB/extern/include' -L'/path_to_MATLAB/bin/maci64' -l'libopencv_core.2.4.13.dylib' -l'libopencv_highgui.2.4.13.dylib' -l'libopencv_imgproc.2.4.13.dylib' -l'libmwlapack.dylib'

---

Troubleshooting Notes:
- The opencv libraries can be difficult for MATLAB to find, even if the linking commands are correct. If it's not finding opencv, you can add the paths to opencv to your PATH, which could help
  - In MATLAB terminal, this would be <code>addpath "/path_to_opencv/stuff"</code>
- The libraries have different names in Windows and macOS, and on top of that, Windows libs end in .lib, while in macOS they end in .dylib
- Libraries can often be in different places in Windows and macOS versions of MATLAB. For instance, the libmwlapack library is at "\path_to_MATLAB\extern\lib\win64\", while in macOS it's at "/path_to_MATLAB/bin/maci64/libmwlapack.dylib"
- The <code>-v</code> flag indicates a verbose arguement, which can be helpful if it doesn't compile correctly. 
- A potential issue with this command that I ran into on macOS is that the libraries weren't linking properly, and I had to explictly define them. Thus, instead of an arguement like <code>-l'libopencv_core.2.4.13.dylib</code>, you would have <code>'/explicit_path/libmwlapack.dylib'</code>
- Another problem on macOS was dependency conflicts with incorrect versions of ffmpeg, which is a dependency. If ffmpeg is installed, but opencv isn't finding it because its not the right version, use brew to install the correct version, then make an alias of that version in '/usr/local/opt/' that links to it, while renaming the existing alias of ffmpeg so that it doesn't get lost.
    - For example, I started out with '/usr/local/opt/ffmpeg/' that linked to '/usr/local/Cellar/ffmpeg/5.0/'. However, it wasn't the right version
    - After installing an older version of ffmpeg and modifying the aliases, I had '/usr/local/opt/ffmpeg/' linked to '/usr/local/Cellar/ffmpeg@4/4.4.1/', and '/usr/local/opt/ffmpeg@5/' linked to '/usr/local/Cellar/ffmpeg/5.0/'


## 4. Examples
*Some high-quality ellipse detection examples run with default parameters and on the same computer with Intel Core i7-7500U 2.7GHz CPU and 8 GB memory*

### 4.1 Detecting all ellipses in the image

---
- The number of detected ellipses: 4; Running time: 0.090s; Resolution: 651 x 436
  <img src="./pics/43_result.jpg" width="73%" height="73%">  
  
---
- The number of detected ellipses: 25; Running time: 0.460s; Resolution: 720 x 435
  <img src="./pics/27_result.jpg" width="73%" height="73%"> 


---
- The number of detected ellipses: 3; Running time: 0.060s; Resolution: 512 x 456
  <img src="./pics/23_result.jpg" width="73%" height="73%"> 


---
- The number of detected ellipses: 8; Running time: 0.110s; Resolution: 752 x 525
  <img src="./pics/666_result.jpg" width="73%" height="73%"> 


### 4.2 Detecting the ellipses with positive polarity  
- The number of detected ellipses: 4; Running time: 0.080s; Resolution: 752 x 525
  <img src="./pics/666_positive.jpg" width="73%" height="73%"> 

### 4.3 Detecting the ellipses with negative polarity
- The number of detected ellipses: 4; Running time: 0.086s; Resolution: 752 x 525
  <img src="./pics/666_negative.jpg" width="73%" height="73%"> 

### 4.4 Detecting the ellipses sharing different polarity  
- The number of detected ellipses: 5; Running time: 0.226s; Resolution: 1000 x 680. ($T_{ac} = 165^{o}$, $T_r = 0.5$)  
  <img src="./pics/different-polarity-detection_all.jpg" width="73%" height="73%"> 


## 5. Successful Application Cases Up to Now
- Car Wheel Hub Recognition
- PCB Inspection
- Object Fingerprinting
- Robot Vision


## 6. Citation
```
@article{lu2019arc,
  title={Arc-Support Line Segments Revisited: An Efficient High-Quality Ellipse Detection},
  author={Lu, Changsheng and Xia, Siyu and Shao, Ming and Fu, Yun},
  journal={IEEE Transactions on Image Processing},
  volume={29},
  pages={768--781},
  year={2020},
  publisher={IEEE}
}
```  

## 7. Our Previous Work  
We also proposed a [circle detection method](https://github.com/AlanLuSun/Circle-detection) in our previous work which could detect circles from image efficiently, precisely and robustly.

