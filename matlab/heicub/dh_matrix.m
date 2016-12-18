function G = dh_matrix(a, d, alpha, theta)
 
% Forward Kinematics : The Denavit-Hartenberg Convention
% http://www.cs.duke.edu/brd/Teaching/Bio/asmb/current/Papers/chap3-forward-kinematics.pdf
 
ctheta = cos(theta);
stheta = sin(theta);
calpha = cos(alpha);
salpha = sin(alpha);
g1 = ctheta;
g2 = -stheta*calpha;
g3 = stheta*salpha;
g4 = g1*a;
g5 = stheta;
g6 = ctheta*calpha;
g7 = -ctheta*salpha;
g8 = g5*a;
g10 = salpha;
g11 = calpha;
g12 = d;
G = [
    g1 g2  g3  g4;
    g5 g6  g7  g8;
    0  g10 g11 g12;
    0  0   0   1
    ];
