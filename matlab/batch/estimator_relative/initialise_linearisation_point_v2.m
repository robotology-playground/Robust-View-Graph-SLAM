function [xs, p] = initialise_linearisation_point_v2(options, p, Pkin)
%xs = initialise_linearisation_point(p, P, options)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% configurations
switch_config;

% print out
disp('   ');
if ~isempty(Pkin)
    if USE_VISION
        disp('Pkin is given, and USE_VISION = 1, kinematics->t and vision->R');
    else
        disp('Pkin is given, and USE_VISION = 0, kinematics->t and kinematics->R');
    end
elseif USE_VISION
    disp('Pkin is empty, and USE_VISION = 1, random->t and vision->R');
else
    error('USE_VISION should be set to 1 and/or Pkin should not be empty');
end

% get kinematics
ncams = length(p);
xc = zeros(6*ncams,1);
tnorm = zeros(ncams,1);
vis = zeros(1,length(p(1).x));
for j = 2 : ncams
    
    flag = 0; % Essential matrix - RANSAC success flag
    status = (p(1).s==1)&(p(j).s==1); % status vector, if track is available
    
    % Using kinematics
    if ~isempty(Pkin)
        %x1 = camera_matrix_to_pose(Pkin{p(1).c});
        c1 = camera_matrix_to_pose(Pkin{1});
        %kin = [Pkin{1};0 0 0 1]\[Pkin{i};0 0 0 1];
        %xc((i-1)*6+(1:3))= kin(1:3,4);
        %xc((i-1)*6+(4:6))= R2w(kin(1:3,1:3))';
        %xi = camera_matrix_to_pose(Pkin{p(i).c});
        cj = camera_matrix_to_pose(Pkin{j});
        pose = transform_to_relative_w(cj, c1);
        xc((j-1)*6+(1:3)) = pose(1:3);
        xc((j-1)*6+(4:6)) = pose(4:6);
        
    % Random initialisation (avoid errors of initialising with zeros)
    else
        if j == 2
            xc((j-1)*6+(1:3)) = [68; 0; 0]/1000;
            if mod(p(j).c,2) % odd camera, even to odd constraint
                xc((j-1)*6+1) = -xc((j-1)*6+1); % reverse baseline sign
            end
        else
            xc((j-1)*6+(1:3))= 1e-3;
            %xc((i-1)*6+(1:3))= rand(3,1)/99;
        end
        xc((j-1)*6+(4:6))= 1e-3;
        %xc((i-1)*6+(4:6))= rand(3,1)*pi/180;
        
    end
    
    % Using vision: Essential matrix with RANSAC
    % Avoid geometries with few matches
    if USE_VISION==1 && sum(status)>options.mincorrnr && (STEREO==0||j>2)
        pwg.edge = [p(1).c, p(j).c];
        use = find(status==1);
        x1 = [p(1).x; p(1).y];
        x2 = [p(j).x; p(j).y];
        if 1
            clf; plot(x1(1,use),x1(2,use),'o');
            hold on; plot(x2(1,use),x2(2,use),'ro');
            title(num2str([p(1).c, p(j).c])); pause;
        end
        pwg = test_Emat_v3(pwg, x1, x2, use, options);
        if sum(pwg.maxinlier)>options.mininlnr % Avoid geometries with few inliers
            flag = 1;
            %xc((j-1)*6+(1:3))= pwg.t*.068;
            xc((j-1)*6+(4:6))= R2w(pwg.r)';
            p(j).s(use) = pwg.maxinlier;
            vis(use) = vis(use)|pwg.maxinlier; % using Essential-RANSAC
        end
    end
    
    % Geometry can not be computed, then random-start
    if flag==0
        vis(status==1) = 1; % using both, Pkin and random initialisation
    end
        
    % Compute the norm of translation
    %tnorm(i) = norm(xc((i-1)*6+(1:3))) ; % 3D displacement vector
    R = w2R(xc((j-1)*6+(4:6)));
    tnorm(j) = norm(R(:,1)'*xc((j-1)*6+(1:3))); % horizontal displacement only
    
    % Print out
    if 1;%j == 2
        fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d\n', ...
            xc((j-1)*6+(1:3))', xc((j-1)*6+(4:6))'*180/pi, sum(p(j).s));
    end
    
end

if USE_VISION==1% remove outliers from the tracks
    for i = 1 : length(p)
        p(i).x = p(i).x(vis==1);
        p(i).y = p(i).y(vis==1);
        p(i).s = p(i).s(vis==1);
    end
    fprintf('%d tracks were removed by RANSAC.\n',sum(vis==0));
end

% assign track visibility based on baseline
npts = length(p(1).s);
imax = zeros(1,npts); % id of camera with largest base-line
%tmax = options.minbase*ones(1, npts) ; % sets the minimum baseline to minbase
tmax = zeros(1,npts); % sets the minimum baseline to zero
tmat = zeros(ncams,npts); % visibility matrix (t-matrix, defined using t)
for i = 1 : ncams
    vis = zeros(1,npts);
    % find a point viewed by a larger baseline than before
    %vis(p(i).s==1) = tnorm(i)*ones(1, sum(p(i).s==1)) > tmax(p(i).s==1);
    vis(p(i).s==1) = ( tnorm(i) > tmax(p(i).s==1) );
    imax(vis==1) = i; % which camera ?
    tmax(vis==1) = tnorm(i); % what baseline ?
    tmat(i,p(i).s==1) = tnorm(i); % visibility is defined by the baseline (not 1/0)
end

% remove points viewed from very short baselines or not visible in any
% camera. This means that imax can not be 0 or 1 (i.e. not seen or only
% visible in the first camera).
use = imax > 1;
for i = 1 : length(p)
    p(i).x = p(i).x(use==1);
    p(i).y = p(i).y(use==1);
    p(i).s = p(i).s(use==1);
end
imax = imax(1, use==1);
fprintf('%d tracks were removed for not being seen by at least %d cameras.\n',sum(use==0),2);

% triangulate points
npts = sum(p(1).s); % new number of points
xf = zeros(npts, 1);
vis = p(1).s;
for i = 1 : npts
    
    x1 = [p(1).x(i);p(1).y(i)];
    %idx = imax(i); % use camera with the largest baseline ?
    %idx = 2 ; % use camera 2 ?
    if p(2).s(i)==1;
        idx = 2;
    else
        idx = imax(i);
    end
    x2 = [p(idx).x(i);p(idx).y(i)];
    
    % consider degeneracy while choosing the second camera frame ?
    %[~,I] = sort(tmat(:,i),1,'descend'); k=1; idx=I(k);
    %x2 = [p(idx).x(i);p(idx).y(i)]; %use camera with the largest baseline ?
    %while isdegenerate(x1, x2)
    %    k = k+1; idx = I(k); %use camera with the (next) largest baseline ?
    %    x2 = [p(idx).x(i); p(idx).y(i)];
    %end
    
    % inverse depth
    xf(i) = test_triangulate_inverse_depth(x1, x2, xc((idx-1)*6+(1:6)));
    
    % remove negative inverse depth ?
    if (xf(i) <= 0); vis(i) = 0; end
    
end

% remove points with negative depth
xf = xf(vis == 1);
for i = 1 : length(p)
    p(i).x = p(i).x(vis==1);
    p(i).y = p(i).y(vis==1);
    p(i).s = p(i).s(vis==1);
end
fprintf('%d tracks were removed for having negative depth.\n',sum(vis==0));

% state and switch vectors
xs = [xc; xf(:)];
%for i = 1 : ncams
%    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d\n', ...
%        xs((i-1)*6+(1:3))', xs((i-1)*6+(4:6))'*180/pi, sum(p(i).s));
%end
%pause

% Check data associations between frames 1 and 2
if 0
    
    status = p(1).s==1 & p(2).s==1;
    im1=get_image_k(options,1);
    im2=get_image_k(options,2);
    clf; imshow([im1 im2]); hold on;
    
    [K,kc]=get_intrinsics(options,1);
    [x,y]=UndistPointFor(p(1).x(status==1),p(1).y(status==1),...
        kc(1), kc(2), kc(3), kc(4), kc(5));
    x=K(1,1)*x+K(1,3);
    y=K(2,2)*y+K(2,3);
    p1=[x;y];
    plot(p1(1,:),p1(2,:),'o');
    
    [K,kc]=get_intrinsics(options,2);
    [x,y]=UndistPointFor(p(2).x(status==1),p(2).y(status==1),...
        kc(1), kc(2), kc(3), kc(4), kc(5));
    x=K(1,1)*x+K(1,3);
    y=K(2,2)*y+K(2,3);
    p2=[x+size(im1,2);y];
    plot(p2(1,:),p2(2,:),'o');
    
    line([p1(1,:);p2(1,:)],[p1(2,:);p2(2,:)],'color',[0 1 0]);
    pause
    
end

% plots
if 0
    r = 1./xf'; % range from inverse range
    X = get_scan_from_range([p(1).x;p(1).y],r);
    clf; plot3(X(1,:),X(3,:),-X(2,:),'+'); axis equal;
    pause;
end
if 0
    im = get_image_k(options,k);
    clf; imshow(im); hold on;
    [K, kc] = get_intrinsics(options,k);
    x1 = [p(1).x;p(1).y];
    kpts = add_lens_distortion(x1,kc);
    kpts = pflat(K*pextend(kpts(1:2,:)));
    plot(kpts(1,:),kpts(2,:),'+') ; axis tight;
    pause;
end
%
%
function x = camera_matrix_to_pose(P)
x = [P(1:3,4); R2w(P(1:3,1:3))'];
%
%
function [p, dp] = get_scan_from_range(p, r)
if size(p,1) < 3; p = pextend(p); end;
rim = sqrt(p(1,:).*p(1,:) + p(2,:).*p(2,:) + p(3,:).*p(3,:));
d = r./rim;
p(1,:) = d.*p(1,:);
p(2,:) = d.*p(2,:);
p(3,:) = d.*p(3,:);
if nargout>1
    dd = -r.*d;
    dp = p;
    dp(1,:) = dd.*p(1,:);
    dp(2,:) = dd.*p(2,:);
    dp(3,:) = dd.*p(3,:);
end
%
%
function im = get_image_k(options, k)
if mod(k,2);
    I = imread(strcat(options.img_folder{1}, options.cam_left.image{(k+1)/2}));
    if size(I,3) == 1 % bayer decoding?
        I = demosaic(I, 'grbg');
    end
    %I = I(5:end-4,5:end-4,:);
else
    I = imread(strcat(options.img_folder{2}, options.cam_right.image{(k+0)/2}));
    if size(I,3) == 1 % bayer decoding?
        I = demosaic(I, 'grbg');
    end
end
if length(size(I)) == 3
    %I = I(:,:,2);
    I = rgb2gray(I);
end
im = I;
%
%
function r = isdegenerate(x1,x2)
r = iscolinear(x1(1,:), x1(2,:), x2(1,:)) || ...
    iscolinear(x1(1,:), x1(2,:), x2(2,:)) || ...
    iscolinear(x1(1,:), x2(1,:), x2(2,:)) || ...
    iscolinear(x1(2,:), x2(1,:), x2(2,:));
%
%
function r = iscolinear(p1, p2, p3, flag)

if nargin == 3   % Assume inhomogeneous coords
    flag = 'inhomog';
end

if ~all(size(p1)==size(p2)) || ~all(size(p1)==size(p3)) || ...
        ~(length(p1)==2 || length(p1)==3)
    error('points must have the same dimension of 2 or 3');
end

% If data is 2D, assume they are 2D inhomogeneous coords. Make them
% homogeneous with scale 1.
if length(p1) == 2
    p1(3) = 1; p2(3) = 1; p3(3) = 1;
end

if flag(1) == 'h'
    % Apply test that allows for homogeneous coords with arbitrary
    % scale.  p1 X p2 generates a normal vector to plane defined by
    % origin, p1 and p2.  If the dot product of this normal with p3
    % is zero then p3 also lies in the plane, hence co-linear.
    r =  abs(dot(cross(p1, p2),p3)) < eps;
else
    % Assume inhomogeneous coords, or homogeneous coords with equal
    % scale.
    r =  norm(cross(p2-p1, p3-p1)) < eps;
end