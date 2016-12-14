% get rotation and its derivatives ('rodrigues', 'exp_skew', 'dcm')
function [R, dRa, dRb, dRc] = get_rotation_matrix(a, method)

% Carl Olsson's notation of rotation:
% R = exp(skew(a)); % good results with the jacobians
% another way:
% R = expm(Ba*a(1)+Bb*a(2)+Bc*a(3)); % bad results with the jacobians
% or
% R = rodrigues(a); % bad results with the jacobians
% notice that: 
%       exp and expm are different
%           exp = element-by-element exponential
%           expm(X) = V*diag(exp(diag(D)))/V; where, [V,D] = EIG(X);
%           expm is used in rodrigues rotations
%       Ba*a(1)+Bb*a(2)+Bc*a(3) = skew(a)
%           where
%               Ba = [0 1 0; -1 0 0; 0 0 0];
%               Bb = [0 0 1; 0 0 0; -1 0 0];
%               Bc = [0 0 0; 0 0 1; 0 -1 0];

if strcmp(method, 'rodrigues');
    % using rodrigues function
    % R = expm(Ba*xc(4)+Bb*xc(5)+Bc*xc(6)); % check update_var.m
    [R, dR] = rodrigues(a);
    dRa = reshape(dR(:, 1), 3, 3);
    dRb = reshape(dR(:, 2), 3, 3);
    dRc = reshape(dR(:, 3), 3, 3);
    
%     [R, dR] = cv.Rodrigues(a); dR = dR';
%     dRa = reshape(dR(:, 1), 3, 3)';
%     dRb = reshape(dR(:, 2), 3, 3)';
%     dRc = reshape(dR(:, 3), 3, 3)';
    
elseif strcmp(method, 'exp_skew');
    % using exp(skew) formula
    % base tangent plane to the rotational diversity
    % rotation derivatives when R = exp(skew(a))
    % R = exp(Ba*xc(4)+Bb*xc(5)+Bc*xc(6));
    R = exp(skew(a));
    Ba = [0 1 0; -1 0 0; 0 0 0];
    Bb = [0 0 1; 0 0 0; -1 0 0];
    Bc = [0 0 0; 0 0 1; 0 -1 0];
    dRa = Ba.*R;
    dRb = Bb.*R;
    dRc = Bc.*R;
elseif strcmp(method, 'dcm');
    % using direction cosine matrix
    R = euler_to_R(a);
    [dRa, dRb, dRc] = euler_to_dR(a);
end