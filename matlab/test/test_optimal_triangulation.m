function X = test_optimal_triangulation(t, R, x1, x2)

if size(x1, 1) < 3;
    x1 = extend_points(x1);
    x2 = extend_points(x2);
end

nbr = size(x1, 2);
X = zeros(4, nbr);

% projection matrices
P1 = eye(3,4);
P2 = [R' -R'*t];

% fundamental matrix
F = ptobi(P1, P2);
F = F/norm(F);
e1 = pflat(cross(F(:,1), F(:,2)));
e2 = pflat(cross(F(1,:)', F(2,:)'));

T1 = eye(3); T1(:,3) = e1;
T2 = eye(3); T2(:,3) = e2;

Ftmp = T1'*F*T2;
[u, s, v] = svd(Ftmp(1:2, 1:2));
T1(1:2, 1:2) = T1(1:2,1:2)*u';
T2(1:2, 1:2) = T2(1:2,1:2)*v;
Ftmp = T1'*F*T2;
lambda = Ftmp(1,1)/Ftmp(2,2);

x1 = T1\x1;
x2 = T2\x2;
P1 = T1\P1;
P2 = T2\P2;

M = [P1, zeros(3,2); P2, zeros(3,2)];
M(3,5) = 1; 
M(6,6) = 1;

for ii = 1:nbr
    
    mx1 = x1(1, ii);
    my1 = x1(2, ii);
    mx2 = x2(1, ii);
    my2 = x2(2, ii);
    
    %syms s;dist=(s*mx1+my1)^2/(s^2+1)+(lambda*mx2-s*my2)^2/(s^2+lambda^2)
    %deriv=diff(dist,s);
    %n=maple('numer',diff(dist,s))
    %maple('coeff',n,'s',6)/2
    %syms s mx1 mx2 my1 my2 lambda
    
    poly=[mx1*my1-mx2*lambda*my2,...
         -lambda^2*my2^2+mx2^2*lambda^2-mx1^2+my1^2,...
         -mx1*my1+2*mx1*my1*lambda^2+mx2*lambda^3*my2-2*mx2*lambda*my2,...
         -2*lambda^2*my2^2-2*mx1^2*lambda^2+2*mx2^2*lambda^2+2*my1^2*lambda^2,...
          mx1*my1*lambda^4-2*mx1*my1*lambda^2+2*mx2*lambda^3*my2-mx2*lambda*my2,...
         -mx1^2*lambda^4-lambda^2*my2^2+mx2^2*lambda^2+my1^2*lambda^4,...
         -mx1*my1*lambda^4+mx2*lambda^3*my2];
    
    r = roots(poly);
    r = r(imag(r) == 0);
    
    [dist, ind] = min((r*mx1+my1).^2./(r.^2+1)+(lambda*mx2-r*my2).^2./(r.^2+lambda^2));
    
    %noise free points
    s = r(ind);
    
    n = [s; 1]/sqrt(s^2+1);
    tmp = [mx1; my1] - [-1; s];
    nx1 = [mx1; my1] - tmp'*n*n;
    
    n = [lambda; -s]/sqrt(s^2 + lambda^2);
    tmp = [mx2; my2] - [s; lambda];
    nx2 = [mx2; my2] - tmp'*n*n;
    
    M(1:2, 5) = nx1;
    M(4:5, 6) = nx2;
    [u,s,v] = svd(M);
    X(:,ii) = pflat(v(1:4,end));
    
    %tmp = null(M);
    %X(:,ii) = pflat(tmp(1:4));
    
end

X = X(1:3, :);


function B=ptobi(P1,P2)
% B=ptobi(m,index)
% Calculates the bilinear tensor from the two camera matrices in MOTION m

for i = 1:3
    for j = 1:3
        B(i,j) = det([P1(1+mod(i,3),:) ; P1(1+mod(i+1,3),:) ; P2(1+mod(j,3),:) ; P2(1+mod(j+1,3),:)]);
    end
end