%function [SOLS,EE] = fivePoint(Q1,Q2)
%
% Copyright Chris Engels 2004
%
% The algorithm follows 
%@Article{         nister-itpam-04,
%  author        = {Nist\'er, D.},
%  journal       = pami,
%  month         = {June},
%  number        = {6},
%  title         = {Problem},
%  pages         = {756-770},
%  volume        = {26},
%  year          = {2004}
%}
%
% This implemenation was written for the paper 
%
%  ARTICLE{stewenius-engels-nister-isprsj-2006,
%  AUTHOR = {H. Stew\'enius and C. Engels and D. Nist\'er},
%  TITLE = {Recent Developments on Direct Relative Orientation},
%  JOURNAL = {ISPRS Journal of Photogrammetry and Remote Sensing},
%  URL = {http://dx.doi.org/10.1016/j.isprsjprs.2006.03.005},
%  VOLUME = {60},
%  ISSUE = {4},
%  PAGES = {284--294},
%  MONTH = JUN,
%  CODE = {http://vis.uky.edu/~stewe/FIVEPOINT},
%  PDF = {http://www.vis.uky.edu/~stewe/publications/stewenius_engels_nister_5pt_isprs.pdf},
%  YEAR = 2006
%}
%
%
% Please refer to this paper if you use this code. 
function [SOLS, EE] = calibrated_fivepoint_non_gb(Q1, Q2)

Q1 = Q1';
Q2 = Q2';

Q = [Q1(:,1).*Q2(:,1) , ...
     Q1(:,2).*Q2(:,1) , ...
     Q1(:,3).*Q2(:,1) , ... 
     Q1(:,1).*Q2(:,2) , ...
     Q1(:,2).*Q2(:,2) , ...
     Q1(:,3).*Q2(:,2) , ...
     Q1(:,1).*Q2(:,3) , ...
     Q1(:,2).*Q2(:,3) , ...
     Q1(:,3).*Q2(:,3) ] ; 

% EE = null(Q); 
[U, S, V] = svd(Q);
EE = V(:,6:9);


A = calibrated_fivepoint_helper( EE );


p=[1 4 2 3 5 11 7 13 6 12 8 14 17 9 15 18 10 16 19 20] ;%rearrange the columns of A
A=A(:,p);
A=rref(A) ;

B(1,1) =-A(6,11);
B(1,2) = A(5,11)-A(6,12);
B(1,3) = A(5,12)-A(6,13);
B(1,4) = A(5,13);
B(1,5) =-A(6,14);
B(1,6) = A(5,14)-A(6,15);
B(1,7) = A(5,15)-A(6,16);
B(1,8) = A(5,16);
B(1,9) =-A(6,17);
B(1,10)= A(5,17)-A(6,18);
B(1,11)= A(5,18)-A(6,19);
B(1,12)= A(5,19)-A(6,20);
B(1,13)= A(5,20);

B(2,1) =-A(8,11);
B(2,2) = A(7,11)-A(8,12);
B(2,3) = A(7,12)-A(8,13);
B(2,4) = A(7,13);
B(2,5) =-A(8,14);
B(2,6) = A(7,14)-A(8,15);
B(2,7) = A(7,15)-A(8,16);
B(2,8) = A(7,16);
B(2,9) =-A(8,17);
B(2,10)= A(7,17)-A(8,18);
B(2,11)= A(7,18)-A(8,19);
B(2,12)= A(7,19)-A(8,20);
B(2,13)= A(7,20);

B(3,1) =-A(10,11);
B(3,2) = A(9,11)-A(10,12);
B(3,3) = A(9,12)-A(10,13);
B(3,4) = A(9,13);
B(3,5) =-A(10,14);
B(3,6) = A(9,14)-A(10,15);
B(3,7) = A(9,15)-A(10,16);
B(3,8) = A(9,16);
B(3,9) =-A(10,17);
B(3,10)= A(9,17)-A(10,18);
B(3,11)= A(9,18)-A(10,19);
B(3,12)= A(9,19)-A(10,20);
B(3,13)= A(9,20);

b11=B(1,1:4)  ;
b12=B(1,5:8)  ;
b13=B(1,9:13) ;

b21=B(2,1:4)  ;
b22=B(2,5:8)  ;
b23=B(2,9:13) ;

b31=B(3,1:4)  ;
b32=B(3,5:8)  ;
b33=B(3,9:13);

n=(conv(conv(b11,b22),b33) - conv(conv(b11,b23),b32)) +...
  (conv(conv(b12,b23),b31) - conv(conv(b12,b21),b33)) +...
  (conv(conv(b13,b21),b32) - conv(conv(b13,b22),b31));

SOLS(:,3)=roots(n);
r=SOLS(:,3);

for i=1:length(r)
    if isreal(r(i))
        
        bt(1,1)=b11(1)*r(i)^3 + b11(2)*r(i)^2 + b11(3)*r(i)^1 + b11(4);
        bt(1,2)=b12(1)*r(i)^3 + b12(2)*r(i)^2 + b12(3)*r(i)^1 + b12(4)  ;
        bt(1,3)=b13(1)*r(i)^4 + b13(2)*r(i)^3 + b13(3)*r(i)^2 + b13(4)*r(i)^1 + b13(5);
        bt(2,1)=b21(1)*r(i)^3 + b21(2)*r(i)^2 + b21(3)*r(i)^1 + b21(4)  ;
        bt(2,2)=b22(1)*r(i)^3 + b22(2)*r(i)^2 + b22(3)*r(i)^1 + b22(4)  ;
        bt(2,3)=b23(1)*r(i)^4 + b23(2)*r(i)^3 + b23(3)*r(i)^2 + b23(4)*r(i)^1 + b23(5);
        bt(3,1)=b31(1)*r(i)^3 + b31(2)*r(i)^2 + b31(3)*r(i)^1 + b31(4)  ;
        bt(3,2)=b32(1)*r(i)^3 + b32(2)*r(i)^2 + b32(3)*r(i)^1 + b32(4)  ;
        bt(3,3)=b33(1)*r(i)^4 + b33(2)*r(i)^3 + b33(3)*r(i)^2 + b33(4)*r(i)^1 + b33(5);

        [U,S,V] = svd( bt);
        xy1 = V(:,end);
        if( xy1(3) ~= 0)
            SOLS(i,1)=xy1(1)./xy1(3);
            SOLS(i,2)=xy1(2)./xy1(3);
        end
    end
end

% %% added by tariq
%
% EE = EE*[SOLS' ; ones(1,10 ) ];
% EE = EE./ ( ones(9,1)*sqrt( sum( EE.^2 ) ) );
%
% %I = find(not(imag( Evec(1,:) )));
% I = not(imag( EE(1,:) ));
% EE = EE(:,I);