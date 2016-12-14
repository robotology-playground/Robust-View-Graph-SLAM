% Computes fundamental matrix from 3 points (Todd's method)
%
%
% Arguments:
%          x1, x2 - Two sets of corresponding 3xN set of homogeneous
%          points.
%
%          x      - If a single argument is supplied it is assumed that it
%                   is in the form x = [x1; x2]
% Returns:
%          F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.


function F = FMatrix_3pts(varargin)

arg  = varargin(:);
x1   = arg{1}(1:3,:)
x2   = arg{1}(4:6,:)
R    = arg{2}
K    = arg{3}
npts = size(x1,2);


% Normalise each set of points so that the origin
% is at centroid and mean distance from origin is sqrt(2).
% normalise2dpts also ensures the scale parameter is 1.
%[x1, T1] = normalise2dpts(x1);
%[x2, T2] = normalise2dpts(x2);
origin = [K(1,2); K(2,3)];
[x1, x2, T1, T2] = Normalization(x1, x2, origin);

A = T2*K*inv(R)*inv(K)*inv(T1)

a11=A(1,1); a12=A(1,2); a13=A(1,3);
a21=A(2,1); a22=A(2,2); a23=A(2,3);
a31=A(3,1); a32=A(3,2); a33=A(3,3);
b = zeros(3,1);
c = zeros(3,1);
d = zeros(3,1);
for i = 1:size(x1,2);
    b(i) = a22*x1(2,i)+a23+a21*x1(1,i)-a33*x2(2,i)-a32*x1(2,i)*x2(2,i)-a31*x1(1,i)*x2(2,i);
    c(i) = -a13-a11*x1(1,i)-a12*x1(2,i)+a33*x2(1,i)+a32*x1(2,i)*x2(1,i)+a31*x1(1,i)*x2(1,i);
    d(i) = -a21*x1(1,i)*x2(1,i)-a22*x1(2,i)*x2(1,i)-a23*x2(1,i)+a12*x1(2,i)*x2(2,i)+a13*x2(2,i)+a11*x1(1,i)*x2(2,i);
end

% Build the constraint matrix
BCD = [b(1) c(1) d(1);
       b(2) c(2) d(2);
       b(3) c(3) d(3)]
opts.UT = true; opts.TRANSA = true;
[e, r] = linsolve(BCD,zeros(3,1), opts)

[U,D,V] = svd(A); % Under MATLAB use the economy decomposition
V(:,end)

ex = [0,   -e(3), e(2);
      e(3), 0,   -e(1);
     -e(2), e(1), 0];
 
F = A*ex;

% Denormalise
F = T2'*F*T1;

% impose the condition of ||F|| = 1
F = F/norm(F,2)

