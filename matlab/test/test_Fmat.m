

% fundamental matrix estimation with MAPSAC.
% Tariq Abuhashim - August 2014, iCub

function [F, cpts1, cpts2] = test_Fmat(mpts1, mpts2, siz, thres, numsmp, T)


% thres : threshold on the norm of the error
% numsmp : number of point samples to be used to calculate F
% T : mapsac threshold

% print text
fprintf('F mapsac .... '); t_start = tic;

% defaults
if nargin < 4;
    thres = .2; % 2 with stereo setting, .2 with monocular
    numsmp = 100; % was 100
    T = .1; % 1 with stereo images, .1 (or 1 with monocular if matching THRESH = 200)
end

error = true;
itr = 0;
while error;
    
    itr = itr + 1;
    
    % the points
    x1 = mpts1(1, :)'; y1 = mpts1(2, :)';
    x2 = mpts2(1, :)'; y2 = mpts2(2, :)';
    
    %[K_left, K_right]=test_param( );
    %mpts1_n = normalise_points(mpts1', K_left);
    %mpts2_n = normalise_points(mpts2', K_right);
    
    %x1 = mpts1_n(1, :)'; y1 = mpts1_n(2, :)';
    %x2 = mpts2_n(1, :)'; y2 = mpts2_n(2, :)';
    
    % fundamental matrix with mapsac
    [f, ~, numin, matches] = torr_mapsac_F(x1, y1, x2, y2, size(x1, 1), 1, numsmp, T);
    F_mapsac = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
    
    % fundamental matrix with ransac
    %[F_mapsac, matches] = ransacfitfundmatrix(mpts1', mpts2', T);
    %numin = size(matches, 2);
    %F = F_mapsac'; f = F(:);
    
    % correct matches.
    x1 = mpts1(1, matches)'; y1 = mpts1(2, matches)';
    x2 = mpts2(1, matches)'; y2 = mpts2(2, matches)';
    [cpts, e] = torr_correctx4F(f, x1, y1, x2, y2, numin, 1);
    cpts = [mpts1(:, matches)' mpts2(:, matches)'];
    
    % check if epipoles are in the image frames
    [e1, e2] = e1_e2_from_F( F_mapsac );
    status = check_epipoles( e1, e2, siz ); % error if status = 1

    % check epipolar geometry status
    if ( norm(e) < thres ) && ~status; error = false; 
    end % epipolar geometry solved, leave loop 
    if itr == numsmp; error = false; cpts = zeros(1, 4);
        fprintf('max number of itr reached, ');
    end % No epipolar geometry found, leave loop
    
end

% outputs
cpts1 = cpts(:, [1 2])'; cpts2 = cpts(:, [3 4])';

if size(cpts, 1) > 7;
    
    % calculate sampson's weight for a set of matches and an F
    x1 = cpts(:, 1); y1 = cpts(:, 2); x2 = cpts(:, 3); y2 = cpts(:, 4);

    % recalculate the fundamental matrix using all the matches without outliers
    fprintf('F all inliers .... '); t_start = tic;
    A(:,1) = x1(:).*x2(:); A(:,2) = y1(:).*x2(:); A(:,3) = x2(:);
    A(:,4) = x1(:).*y2(:); A(:,5) = y1(:).*y2(:); A(:,6) = y2(:);
    A(:,7) = x1(:);        A(:,8) = y1(:);        A(:,9) = 1;
    [U,D,V] = svd(A);
    f = V(:,length(V));
    F = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
    %error = D(length(V),length(V));
    
    % Modify F to be rank-2
    [FU,FD,FV] = svd(F);
    FDnew = FD;
    FDnew(3,3) = 0;
    F_linear = FU*FDnew*FV';

    % which fundamental matrix?
    %F = F_mapsac;
    F = F_linear;
    
else
    
    F = F_mapsac;
    
end

% print text
fprintf([num2str(itr),' itrs, ',num2str(numin),' inliers, ', ...
    num2str(norm(e)),' error, ']);
t_end = toc(t_start); fprintf([num2str(t_end),'s\n']);