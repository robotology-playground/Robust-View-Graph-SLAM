

% fundamental matrix estimation with MAPSAC.
% Tariq Abuhashim - August 2014, iCub

function [E,matches] = test_Fmat_2(p1,p2,K1,K2,opt)


thres = .2; % 2 with stereo setting, .2 with monocular
numsmp = 100; % was 100
T = .1; % 1 with stereo images, .1 (or 1 with monocular if matching THRESH = 200)
error = true;
itr = 0;

while error;
    
    itr = itr + 1;
    
    % the points
    x1=p1(1, :)';y1=p1(2, :)';
    x2=p2(1, :)';y2=p2(2, :)';

    % fundamental matrix with mapsac
    [f, e, ~, matches] = torr_mapsac_F(x1, y1, x2, y2, size(x1,1),1,numsmp,T);
    F_mapsac = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
    
    % fundamental matrix with ransac
    %[F_mapsac, matches] = ransacfitfundmatrix(mpts1', mpts2', T);
    %numin = size(matches, 2);
    %F = F_mapsac'; f = F(:);
    
    % get the matches
    p1=p1(:,matches);p2=p2(:,matches);
    
    % check if epipoles are in the image frames
    [e1,e2]=e1_e2_from_F(F_mapsac);
    status=check_epipoles(e1,e2,opt.imgsize); % error if status = 1

    % check epipolar geometry status
    if (norm(e)<thres)&&~status;error=false; 
    end % epipolar geometry solved, leave loop 
    if itr==numsmp;error=false;
        fprintf('max number of itr reached, ');
    end % No epipolar geometry found, leave loop
    
end


if 0;%size(p1,2)>7;
    
    % calculate sampson's weight for a set of matches and an F
    x1=p1(1, :)';y1=p1(2, :)';
    x2=p2(1, :)';y2=p2(2, :)';
    % recalculate the fundamental matrix using all the matches without outliers
    %fprintf('F all inliers .... '); t_start = tic;
    A(:,1)=x1(:).*x2(:);A(:,2)=y1(:).*x2(:);A(:,3)=x2(:);
    A(:,4)=x1(:).*y2(:);A(:,5)=y1(:).*y2(:);A(:,6)=y2(:);
    A(:,7)=x1(:);       A(:,8)=y1(:);       A(:,9)=1;
    [U,D,V]=svd(A);
    f=V(:,length(V));
    F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
    %error = D(length(V),length(V));
    
    % Modify F to be rank-2
    [FU,FD,FV]=svd(F);
    FDnew=FD;
    FDnew(3,3)=0;
    F_linear=FU*FDnew*FV';

    % which fundamental matrix?
    %F=F_mapsac;
    F=F_linear;
    
else
    
    F = F_mapsac;
    
end

% essential matrix
E=K1'*F*K2;
E=E/norm(E);
[U,S,V]=svd(E);
if abs(S(3,3))>0.001
    error('F must be rank 2 to self calibrate');
end
%S(3,3)=0;
%S(1,1)=S(2,2);
%E=U*S*V';
