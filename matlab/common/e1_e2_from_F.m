


% compute the two epipoles
% Tariq Abuhashim - August 2014, iCub

function [e1, e2] = e1_e2_from_F(F)


% the epipole e1 is the projection of optical centre C2 onto image 1
% the epipole e2 is the projection of optical centre C1 onto image 2

% calculate epipoles
[U,D,V] = svd(F);
e1 = V(:,3); % e1=null(F); % null uses svd too
e2 = U(:,3); % e2=null(F');

% Put epipoles in front of camera
if e1 < 0; e1 = -e1; end
if e2 < 0; e2 = -e2; end