function F = F_from_P1_P2(P1, P2)

% calculates the fundamental matrix from camera matrices P1 and P2
%
% notes:
% notice that if using normalised image coordinates:
%   x_n = K\x
%   then E = F. In this case, the following two functions should produce the
%   same results
%       F = F_from_P1_P2(P1, P2)
%       E = E_from_R_t(R, t)
%
% In order to calculate the fundamental matrix from the essential matrix, 
% use this function
%       F = F_from_E (E, K1, K2)
%   it produces more accurate results than calculating F from normalised
%   camera matrices P1 and P2
%
% Tariq Abuhashim

if ~all(size(P1) == [3 4]) || ~all(size(P2) == [3 4])
    error('Camera matrices must be 3x4');
end

C1 = null(P1);  % Camera centre 1 is the null space of P1
e2 = P2*C1;     % epipole in camera 2

F = skew_r(e2)*P2*pinv(P1);