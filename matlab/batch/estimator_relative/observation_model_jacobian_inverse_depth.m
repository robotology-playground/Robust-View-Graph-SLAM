function H = observation_model_jacobian_inverse_depth(x, p1)
%H = observation_model_jacobian_inverse_depth(x, p1)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% % Code to test the against the numerical jacobian
% addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/matlab/utilities');
% % test data
% p1 = rand(2,1); x = rand(7,1);
% H = observation_model_jacobian_inverse_depth(x, p1);
% J = numerical_jacobian_i(@observation_model_inverse_depth, [], 1, [], x, p1);
% numerical = J
% analytical = full(H)
% error = analytical - numerical

% get rotation and its derivatives ('rodrigues', 'exp_skew', 'dcm')
[RR, dR1, dR2, dR3] = get_rotation_matrix(x(4:6), 'rodrigues');

% Allocate storage for triplet form
N = size(p1,2);
Hii = zeros(2*N, 1); % 2 diagonal terms for each point
Hij = zeros(12*N, 3); % 12 motion terms for each point
k = 1:12;
for i2 = 1:N
    i1 = 1:6;
    ii = [i1 i2+6]';
    ih = getindex_point_triplet(i2);
    hi = compute_observation_model_derivatives(x(ii), p1(:,i2), RR, dR1, dR2, dR3);
    Hii(ih) = hi(1:2);
    % Store block-diagonal terms in triplet form (terms are unique)
    ii2 = getindex2(i2);
    Hij(k,:) = off_diagonal_triplet_form(hi(3:14), ii2);
    k = k + 12;
end
[ii, jj] = compute_block_diagonal_triplet_indices(N);
Hii = sparse(ii, jj, Hii, 2*N, N+6);
Hij = sparse(Hij(:,1), Hij(:,2), Hij(:,3), 2*N, N+6);
H = Hii + Hij;

function idx = getindex2(i)
idx = [2*(i-1) + 1; 2*(i-1) + 2];
idx = idx(:);
function idx = getindex_point_triplet(i)
idx = (1:2) + (i-1)*2;
function Hij = off_diagonal_triplet_form(Hij, i2)
% Get off-diagonal block of Y represented in triplet form
rows = [ones(1,6); 2*ones(1,6)]; rows = rows(:);
cols = [1:6;1:6]; cols = cols(:);
Hij = [rows + i2(1) - 1, cols, Hij(:)];
function [i, j] = compute_block_diagonal_triplet_indices(N)
% Compute indices of the triplet-form block-diagonal terms
i = zeros(2, N);
i(:) = 1:2*N;
j = [1:N; 1:N]+6; % +6 because of 1 camera has 6 parameters
j = j(:);