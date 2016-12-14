% get rotation and its derivatives ('rodrigues', 'exp_skew', 'dcm')
function a = get_rotation_angles(R, method)

if strcmp(method, 'rodrigues');
    % using rodrigues function
    % R = expm(Ba*xc(4)+Bb*xc(5)+Bc*xc(6)); % check update_var.m
    a = rodrigues(R);
elseif strcmp(method, 'exp_skew');
    % using exp(skew) formula
    % base tangent plane to the rotational diversity
    % rotation derivatives when R = exp(skew(a))
    % R = exp(Ba*xc(4)+Bb*xc(5)+Bc*xc(6));
    w = log(R);
    a = [w(1,2); w(1,3); w(2,3)];
elseif strcmp(method, 'dcm');
    % using direction cosine matrix
    a = R_to_euler(R);
end