function [y, Y] = generate_joint_info_matrix(Cinfo, N)
% Form the per-constraint information into a joint {y,Y} where Y is
% sparse. N is the number of poses.

M = length(Cinfo);

% Generate empty information-form
if M == 0
    y = zeros(6*N, 1);
    Y = sparse(6*N, 6*N); 
    return
end

% Allocate storage for triplet form
y = zeros(6*N, 1);   % 3 states for each pose
Yii = zeros(36*N, 1); % 9 diagonal terms for each pose
Yij = zeros(36*M, 3); % 9 off-diagonal terms for each constraint

% Store info in triplet form
k = 1:36;
for i = 1:M
    edge = Cinfo(i).edge;
    i1 = get_state_index(edge(1));
    i2 = get_state_index(edge(2));
    
    % Store info-vector
    ii = [i1 i2];
    y(ii) = y(ii) + Cinfo(i).y;

    % Store block-diagonal terms in triplet form (terms are additive)
    Yi = Cinfo(i).Y;
    iY1 = getindex_triplet(edge(1));
    iY2 = getindex_triplet(edge(2));
    Yii(iY1) = Yii(iY1) + make_col(Yi(1:6,1:6)); % new block-diagonal terms are additive
    Yii(iY2) = Yii(iY2) + make_col(Yi(7:12,7:12));   
    
    % Store off-diagonal terms in triplet form (terms are unique)
    Yij(k,:) = off_diagonal_triplet_form(Yi, i1, i2);
    k = k + 36;
end

% Convert triplet form to sparse col-compressed
[ii, jj] = compute_block_diagonal_triplet_indices(N);
Yii = sparse(ii, jj, Yii, 6*N, 6*N); 
Yij = sparse(Yij(:,1), Yij(:,2), Yij(:,3), 6*N, 6*N);
Y = Yii + Yij + Yij';

%
%

function c = make_col(m)
% Convert matrix into a column vector
c = m(:);

function idx = getindex_triplet(i)
% Get indices into Yii in triplet form
idx = (1:6*6) + (i-1)*6*6;

function Yij = off_diagonal_triplet_form(Y, i1, i2)
% Get off-diagonal block of Y represented in triplet form
rows = repcol([1;2;3;4;5;6], 6); rows = rows(:);
cols = reprow([1,2,3,4,5,6], 6); cols = cols(:);
Yij = Y(1:6, 6+1:6*2);
Yij = [rows + i1(1) - 1, cols + i2(1) - 1, Yij(:)];

function [ii, jj] = compute_block_diagonal_triplet_indices(N)
% Compute indices of the triplet-form block-diagonal terms
ii=zeros(6,N); 
ii(:) = 1:6*N; 
ii = make_col(repmat(ii, 6, 1));
jj = make_col(reprow(1:6*N, 6)); 
