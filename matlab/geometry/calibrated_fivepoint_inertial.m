function Evec = calibrated_fivepoint_inertial(imu)


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

nav_state_prev = eps*ones(9, 1);
nav_state = ins_inertial_mech(nav_state_prev, imu, sim_dt);
t = nav_state(1:3);

%I = find(not(imag(Evec(1,:))));
I = not(imag( Evec(1,:)));
Evec = Evec(:,I);
