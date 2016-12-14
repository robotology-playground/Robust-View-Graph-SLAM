function A = skew_r(a)

% calculates the reversed skew symetric matrix A of elements vector a
% reversed skew matrix (for epipolar geometry stuff)
A = [ 0     -a(3) a(2); 
     a(3)  0     -a(1); 
     -a(2) a(1)   0];