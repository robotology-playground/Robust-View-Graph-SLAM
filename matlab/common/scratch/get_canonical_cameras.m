function [P,Pprime] = get_canonical_cameras( F )

% Get the "canonical" cameras for given fundamental matrix
% according to Hartley and Zisserman (2004), p256, Result 9.14

% But ensure that the left 3x3 submatrix of Pprime is nonsingular
% using Result 9.15, that the general form is
% [ skew( e1 ) * F + e1 * v', k * e1 ] where v is an arbitrary
% 3-vector and k is an arbitrary scalar

P = [ 1 0 0 0;
    0 1 0 0;
    0 0 1 0 ];

e1 = null( F' ); %e1=e1/e1(3);

e1_ = [   0   -e1(3)  e1(2);
       e1(3)   0     -e1(1);
      -e1(2)   e1(1)   0   ];
% e1_ = skew(e1);

M = e1_ * F + e1 * [1 1 1];
Pprime = [ M, e1 ];


% Test that camera matrices Pprime and P are consistent with
% fundamental matrix F
% Meaning  (Pprime*X)' * F * (P*X) = 0,  for all X in 3space
% Get the epipole in camera 1 for camera 2
C2=null(P);
eprime=Pprime*C2;

% Construct F from Pprime, P, and eprime
eprime_ = [   0   -eprime(3)  eprime(2);
       eprime(3)   0     -eprime(1);
      -eprime(2)   eprime(1)   0   ];
% eprime_=skew(eprime);

Fhat = eprime_ * Pprime * pinv(P);
% Check that it's close to F
alpha = Fhat(:)\F(:);
if norm( alpha*Fhat-F ) > 1e-10
    fprintf( 1, 'Warning: supplied camera matrices are inconsistent with F\n' );
end

