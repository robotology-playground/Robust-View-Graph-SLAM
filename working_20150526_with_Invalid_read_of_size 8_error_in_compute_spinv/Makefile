#Manual:
#g++ -D_GNU_SOURCE -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/UFconfig -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/MATLAB -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/AMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/COLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CCOLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/metis-4.0/Lib  -I"/usr/local/MATLAB/R2016a/extern/include" -I"/usr/local/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG main.cpp PwgOptimiser.cpp -o hello


#Makefile:
all: PwgOptimiser

PwgOptimiser: main.o PwgOptimiser.o RecoverMoments.o
	g++ main.o PwgOptimiser.o RecoverMoments.o -o PwgOptimiser -lcholmod -llapack -lblas

main.o: main.cpp
	g++ -c -D_GNU_SOURCE -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/UFconfig -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/MATLAB -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/AMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/COLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CCOLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/metis-4.0/Lib  -I"/usr/local/MATLAB/R2016a/extern/include" -I"/usr/local/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/tabuhashim/Dropbox/code/matlab/cpp/src/main.cpp

PwgOptimiser.o: PwgOptimiser.cpp
	g++ -c -D_GNU_SOURCE -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/UFconfig -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/MATLAB -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/AMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/COLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CCOLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/metis-4.0/Lib  -I"/usr/local/MATLAB/R2016a/extern/include" -I"/usr/local/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/tabuhashim/Dropbox/code/matlab/cpp/src/PwgOptimiser.cpp

RecoverMoments.o: RecoverMoments.cpp
	g++ -c -D_GNU_SOURCE -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/UFconfig -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/MATLAB -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/AMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/COLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CCOLAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CAMD/Include -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/metis-4.0/Lib  -I"/usr/local/MATLAB/R2016a/extern/include" -I"/usr/local/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/tabuhashim/Dropbox/code/matlab/cpp/src/RecoverMoments.cpp

clean:
	rm -r *o PwgOptimiser
