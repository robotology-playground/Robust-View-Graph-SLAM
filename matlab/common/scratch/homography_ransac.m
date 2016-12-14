%
% Calculates the projective transformation between matches keypoints pts1
% and pts2.
% The function outputs, the homography and inliers using ransac
%
% Tariq Abuhashum
% started: 21 August 2014
%
% iCub - Koroibot
%

function [H,inliers]=homography_ransac(pts1,pts2)

% convert to homogeneous coordinates
X1=pts1(1:2,:); X1(3,:)=1;
X2=pts2(1:2,:); X2(3,:)=1;

% ransac with homography model
numMatches=size(pts1,2) ;
for t=1:100
    % estimate homograpyh
    subset=vl_colsubset(1:numMatches,4);
    A=[];
    for i=subset
        A=cat(1,A,kron(X1(:,i)',vl_hat(X2(:,i)))) ;
    end
    [U,S,V]=svd(A) ;
    H{t}=reshape(V(:,9),3,3) ;
    % score homography
    X2_=H{t} * X1 ;
    du=X2_(1,:)./X2_(3,:)-X2(1,:)./X2(3,:);
    dv=X2_(2,:)./X2_(3,:)-X2(2,:)./X2(3,:);
    inliers{t}=(du.*du+dv.*dv)<6*6;
    score(t)=sum(inliers{t});
end

[score,best]=max(score);
H=H{best};
inliers=inliers{best};

% optimal refinement
if exist('fminsearch','file')==2
    H=H/H(3,3) ;
    opts=optimset('Display','none','TolFun',1e-8,'TolX',1e-8);
    H(1:8)=fminsearch(@homography_residual,H(1:8)',opts);
else
    warning('Refinement disabled as fminsearch was not found.');
end

