function Evec = calibrated_fivepoint(Q1, Q2)


%
%  notes:
%  1. reference:
%
%       H. Stew\'enius and C. Engels and D. Nist\'er, 2006
%       Recent Developments on Direct Relative Orientation,
%       http://vis.uky.edu/~stewe/FIVEPOINT,
%       http://www.vis.uky.edu/~stewe/publications/stewenius_engels_nister_5pt_isprs.pdf,
%
%       Henrik Stewenius, 2005
%       Grobner Basis Methods for Minimal Problems in Computer Vision
%       PhD Thesis, Lund University, 2005
%       http://www.maths.lth.se/matematiklth/personal/stewe/THESIS/ 
%
% 2. if this implementation is too slow, please see: 
%
%       Nist\'er, D., 2004
%       An Efficient Solution to the Five-Point Relative Pose
%
% 3. code to veryfy that it works: 
%       Q1 = rand(3,5);
%       Q2 = rand(3,5);
%       Evec   = calibrated_fivepoint( Q1,Q2);
%       for i=1:size(Evec,2)
%           E = reshape(Evec(:,i),3,3);
%           % Check determinant constraint! 
%           det( E)
%           % Check trace constraint
%           2 *E*transpose(E)*E -trace( E*transpose(E))*E
%           % Check reprojection errors
%           diag( Q1'*E*Q2)
%       end
%
% 4. due to varying standards of Q1 and Q2, it is very possible that you 
%    get essential matrices which are the transpose of what your expected.
%

Q1 = Q1';
Q2 = Q2';

Q = [Q1(:,1).*Q2(:,1),...
     Q1(:,2).*Q2(:,1),...
     Q1(:,3).*Q2(:,1),... 
     Q1(:,1).*Q2(:,2),...
     Q1(:,2).*Q2(:,2),...
     Q1(:,3).*Q2(:,2),...
     Q1(:,1).*Q2(:,3),...
     Q1(:,2).*Q2(:,3),...
     Q1(:,3).*Q2(:,3)]; % (5x9) x 1, where 5 points are being used

[~,~,V] = svd(Q, 0); % produces the "economy size" decomposition. 
% If Q is m-by-n with m > n, then only the first n columns of U 
% are computed and S is n-by-n.
% For m <= n, svd(X,0) is equivalent to svd(X).
EE = V(:,6:9); % 10x4
   
M = calibrated_fivepoint_helper(EE);

B = M(:,1:10)\M(:,11:20);
A = zeros(10);
A(1:6,:) = -B([1 2 3 5 6 8], :);
A(7,1) = 1;
A(8,2) = 1;
A(9,4) = 1;
A(10,7) = 1;
[V, ~] = eig(A);
SOLS = V(7:9,:)./(ones(3,1)*V(10,:));

Evec = EE*[SOLS ; ones(1, 10 )]; 
Evec = Evec./ ( ones(9,1)*sqrt( sum( Evec.^2 ) ) );

%I = find(not(imag(Evec(1,:))));
I = not(imag( Evec(1,:)));
Evec = Evec(:,I);
