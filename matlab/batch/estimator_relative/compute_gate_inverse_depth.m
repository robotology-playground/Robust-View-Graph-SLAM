function gate = compute_gate_inverse_depth(x, P, C, xs, ncams)
%gate = compute_gate_inverse_depth(x, P, C, xs)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% Configurations
switch_config;

gate = Inf(1, length(C));
% 36 pose + 1 inverse depth + 6 pose-idepth correlations
%Pii = zeros(43,1);
%ii = zeros(43,1);
%jj = zeros(43,1);
% diagonals (36 pose + 1 inverse depth)
%i = repmat((1:6)',1,6); ii(1:37,1) = [i(:); 7]; % rows
%j = repmat(1:6,6,1); jj(1:37,1) = [j(:); 7]; % cols
% off-diagonals (6 pose-inverse depth correlations higher triangle terms)
%ii(38:43,1) = (1:6)'; % rows
%jj(38:43,1) = (7*ones(1,6))'; % cols
for i = 1:length(C);
    % get indeces
    i1 = (C(i).cam-1)*6 + (1:6)'; % pose index
    i2 = C(i).kpt + ncams*6; % point (also inverse depth) index
    idx = [i1; i2];
    % get marginal covariance
    Pi = P(idx,idx);
    %full(Pi)
    %Pii(1:36,1) = make_col(P(i1,i1));
    %Pii(37,1) = P(i2,i2);
    %Pii(38:43,1) = P(i1,i2);
    %Pi = sparse(ii(1:37),jj(1:37),Pii(1:37),7,7);
    %Pij = sparse(ii(38:43),jj(38:43),Pii(38:43),7,7);
    %Pi = Pi + Pij + Pij';
    % compute the innovations
    zs = observation_model_inverse_depth(xs(idx), C(i).p1);
    Hs = observation_model_jacobian_inverse_depth(xs(idx), C(i).p1);
    dx = x(idx) - xs(idx);
    dx(5) = pi_to_pi(dx(5)); % heading rounds to -pi:pi
    zpred = zs + Hs*dx;
    v = C(i).z - zpred;
    % compute the gate
    S = Hs*Pi*Hs' + sparse(C(i).R);
    warning off MATLAB:nearlySingularMatrix
    gate(i) = v'*(S\v);
    %gate(i) = v'*(inv(S)*v);
    if isnan(gate(i))||~isfinite(gate(i))
        error('rrrrrrrrrrrrr');
    end
    warning on MATLAB:nearlySingularMatrix
    %dx
    %full(Pi)
    %full(Hs)
    %v
    %full(S)
end