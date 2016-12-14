function A = skew(a)

% calculates the skew symetric matrix A of elements vector a

% navigation
A = [ 0     a(1) a(2); 
     -a(1)  0    a(3); 
     -a(2) -a(3) 0];