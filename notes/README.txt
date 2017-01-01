COMPILE INSTRUCTIONS
====================

Change path to libraries and include files accordingly.
lesson5.cpp: generates a Macbeth chart image using spectral data
imageviewer.cpp: a simple PPM image viewer with exposure control

c++ -o lesson5 lesson5.cpp -O3 -Wall
c++ -o imageviewer -L$PATHTOGLFW/glfw-2.7.7/lib/x11 -lglfw -I$PATHTOGLFW/glfw-2.7.7/include -O3 -Wall imageviewer.cpp -I$PATHTOCG/Cg/2.0/linux.centos4.x86_64/include/ -L$PATHTOCG/Cg/2.0/linux.centos4.x86_64/lib64 -lCgGL -lCg
