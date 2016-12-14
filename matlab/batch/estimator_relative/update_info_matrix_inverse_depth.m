function [y, Y] = update_info_matrix_inverse_depth(C, npts, ncams)
%[y, Y] = update_info_matrix_inverse_depth(C, N, ncams)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

M = length(C);
if M == 0;
    y = zeros(npts + 6*ncams, 1);
    Y = sparse(npts + 6*ncams, npts + 6*ncams);
    return;
end

% allocate storage for triple form
y = zeros(6*ncams+npts, 1); % 3 states for each point + 6 for each pose
Yii = zeros(36*ncams+npts, 1); % 9 diaginal terms for each point + 36 for each pose
Yij = zeros(6*M, 3); % 1x6 and 6x1 blocks,

% store information in triplet format
k = 1:6;
for i = 1 : M
    % Store info-vector
    i1 = (C(i).cam-1)*6 + (1:6)'; % pose indeces
    i2 = C(i).kpt + 6*ncams; % inverse depth index
    ii = [i1; i2]; % add pose indeces
    y(ii) = y(ii) + C(i).y;
    % Store block-diagonal terms in triplet form (terms are additive)
    iY1 = getindex_pose_triplet(C(i).cam); % pose indeces
    iY2 = C(i).kpt + 36*ncams; % inverse depth index
    Yi = C(i).Y;
    Yii(iY1) = Yii(iY1) + make_col(Yi(1:6,1:6));
    Yii(iY2) = Yii(iY2) + Yi(7,7);
    % Store block-diagonal terms in triplet form (terms are unique)
    Yij(k,:) = off_diagonal_triplet_form(Yi, i1, i2);
    k = k + 6;
end

% Convert triplet form to sparse col-compressed
[ii, jj] = compute_block_diagonal_triplet_indices(ncams,npts);
Yii = sparse(ii, jj, Yii, 6*ncams+npts, 6*ncams+npts);
Yij = sparse(Yij(:,1), Yij(:,2), Yij(:,3), 6*ncams+npts, 6*ncams+npts);
Y = Yii + Yij + Yij';

function c = make_col(m)
c = m(:);
function idx = getindex_pose_triplet(i)
idx = (i-1)*36+(1:36);
function [i, j] = compute_block_diagonal_triplet_indices(ncams,N)
% Compute indices of the triplet-form block-diagonal terms
% poses
i=[];%(1:6)'; % first camera is diagonal
j=[];%(1:6)'; % first camera is diagonal
for cam = 1:ncams;
    i1 = zeros(6,1);
    i1(:) = ((cam-1)*6)+(1:6);
    i1 = repmat(i1,6,1); i1 = i1(:);
    j1 = reprow(((cam-1)*6)+(1:6),6); j1 = j1(:);
    i = [i; i1];
    j = [j; j1];
end
% 3D points
i2 = zeros(N,1);
i2(:) = 6*ncams+(1:N);
i2 = repmat(i2, 1, 1);
i2 = i2(:);
j2 = reprow(6*ncams+(1:N), 1);
j2 = j2(:);
% full vector
i = [i; i2];
j = [j; j2];
function Yij = off_diagonal_triplet_form(Y, i1, i2)
% Get off-diagonal block of Y represented in triplet form
rows = [1 1 1 1 1 1]; rows = rows(:);
cols = [1 2 3 4 5 6]; cols = cols(:);
Yij = Y(7,1:6);
Yij = [rows + i2(1) - 1, cols + i1(1) - 1, Yij(:)];