function [xs, p] = initialise_linearisation_point(p, Pkin, options, k)
%xs = initialise_linearisation_point(p, P, options)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% configurations
switch_config ;

% print out
disp('   ');
if ~isempty(Pkin)
    if USE_VISION
        disp('Pkin provided, and USE_VISION = 1, kinematics->t and vision->R');
    else
        disp('Pkin provided, and USE_VISION = 0, kinematics->t and kinematics->R');
    end
elseif USE_VISION
    disp('Pkin is empty, and USE_VISION = 1, vision->t and vision->R');
else
    error('USE_VISION should be set to 1 and/or Pkin should not be empty');
end

% get kinematics
ncams = size(p, 1) ;
xc = zeros(6*ncams, 1) ;
tnorm = zeros(ncams, 1) ;
x1 = camera_matrix_to_pose(Pkin{1}) ;
for i = 2 : ncams ;
    
    % using kinematics
    if ~isempty(Pkin)
        %kin = [Pkin{1};0 0 0 1]\[Pkin{i};0 0 0 1];
        %xc((i-1)*6+(1:3))= kin(1:3,4);
        %xc((i-1)*6+(4:6))= R2w(kin(1:3,1:3))';
        xi = camera_matrix_to_pose(Pkin{i}) ;
        pose = transform_to_relative_w(xi, x1) ;
        xc((i-1)*6+(1:3)) = pose(1:3) ;
        xc((i-1)*6+(4:6)) = pose(4:6) ;
    end
    
    % using vision
    if USE_VISION
        
        % reset state to a non-zero (in case no tracks are available)
        % a better solution is to attempt to remove the image before
        % computing the linearisation point .........
        if isempty(Pkin)
            xc((i-1)*6+(1:3))= 1e-4;
        end
        xc((i-1)*6+(1:3))= 1e-4;
        
        % using vision
        status = (p{1}(3,:)==1)&(p{i}(3,:)==1) ; % status vector, if track is available
        if sum(status) > options.mincorrnr % avoid geometries with few matches
            pwg.edge = [1 i]+k-1 ;
            use = find(status) ;
            pwg = test_Emat_v3(pwg, p{1}(1:2,:), p{i}(1:2,:), use, options) ;
            if isempty(Pkin)
                %xc((i-1)*6+(1:3))= pwg.t;%/1000 ;
                xc((i-1)*6+(1:3))= rand(3,1)/100 ;
            end
            xc((i-1)*6+(4:6))= R2w(pwg.r)' ;
        else                             % geometry can not be computed
            if isempty(Pkin)
                xc((i-1)*6+(1:3))= rand(3,1)/100 ;
            end
            xc((i-1)*6+(4:6))= rand(3,1)*pi/180 ;
        end
        if i==2 && isempty(Pkin) % force model second camera displacement
            xc((i-1)*6+(1:3)) = [68; 0; 0]/1000 ;
        end
        
    end
    
    % compute the norm of translation
    %tnorm(i) = norm(xc((i-1)*6+(1:3))) ; % 3D displacement vector
    R = w2R(pose(4:6)) ;
    tnorm(i) = norm(R(:,1)'*xc((i-1)*6+(1:3))) ; % horizontal displacement only
    
    % print out
    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d\n', ...
        xc((i-1)*6+(1:3))', xc((i-1)*6+(4:6))'*180/pi, sum(p{i}(3,:)));
    
end

% assign track visibility based on baseline
npts = size(p{1}, 2) ;
imax = zeros(1, npts) ;
%tmax = options.minbase*ones(1, npts) ; % sets the minimum baseline to minbase
tmax = zeros(1, npts) ; % sets the minimum baseline to zero
tmat = zeros(ncams, npts) ; % visibility matrix (t-matrix, defined using t)
for i = 1 : ncams
    vis = zeros(1, npts) ;
    % find a point viewed by a larger baseline than before
    %vis(1, p{i}(3,:)==1) = tnorm(i)*ones(1, sum(p{i}(3,:)==1)) > tmax(1, p{i}(3,:)==1);
    vis(1, p{i}(3,:)==1) = (tmax(1, p{i}(3,:)==1) < tnorm(i)) ; 
    imax(vis==1) = i ; % which camera ?
    tmax(vis==1) = tnorm(i) ; % what baseline ?
    tmat(i,p{i}(3,:)==1) = tnorm(i) ; % visibility is defined by the baseline (not 1/0)
end

% remove points viewed from very short baselines or not visible in any
% camera (this means that imax can not be 0 or 1)
use = imax > 1 ;
for i = 1 : size(p,1)
    p{i} = p{i}(:, use == 1);
end
imax = imax(1, use == 1) ;
%tmax = tmax(1, use == 1);

% triangulate points
npts = size(p{1}, 2) ; % new number of points
xf = zeros(npts, 1) ;
%sigmas = zeros(npts, 1) ;
vis = zeros(1, npts) ;
%fx = options.K1(1, 1) ;
for i = 1 : npts
    
    x1 = p{1}(1:2, i) ;
    idx = imax(i) ; % use camera with the largest baseline?
    %idx = 2 ; % use camera 2
    x2 = p{idx}(1:2, i);
    
    %[~,I]=sort(tmat(:,i),1,'descend'); k=1; idx=I(k);
    %x2=p{idx}(1:2,i); %use camera with the largest baseline?
    %while isdegenerate(x1,x2)
    %    k=k+1; idx=I(k);
    %    x2=p{idx}(1:2,i); %use camera with the (next) largest baseline?
    %end
    
    % inverse depth
    xf(i) = test_triangulate_inverse_depth(x1, x2, xc((idx-1)*6+(1:6))) ;
    
    % inverse depth uncertainty
    %H = test_triangulate_jacobian_inverse_depth(p{1}(1:2,i), p{imax(i)}(1:2,i), xc((imax(i)-1)*6+(1:6)));
    %sigmas(i) = H*diag(([[0,0]/fx,[3,3]/fx,5.*[1,1,1]/1000,5.*[1,1,1]*pi/180].^2))*H';
    
    % remove negative inverse depth
    %if xf(i) > 0;
        vis(i) = 1;
    %end
    
end

% remove points with negative depth
xf = xf(vis == 1);
%sigmas = sigmas(vis == 1);
for i = 1 : size(p,1)
    p{i} = p{i}(:, vis == 1) ;
end

% plots
if 0
    r = 1./xf'; % range from inverse range
    X = get_scan_from_range(p{1}(1:2,:), r) ;
    clf ; plot3(X(1,:), X(3,:), -X(2,:),'+') ; axis equal ;
    pause ;
end
if 0
    im = get_bundle_images(options, k) ;
    clf; imshow(im{1}) ; hold on ;
    [K, kc] = get_intrinsics(options, k) ;
    kpts = add_lens_distortion(p{1}(1:2,:), kc) ;
    kpts = pflat(K*pextend(kpts(1:2, :))) ;
    plot(kpts(1,:), kpts(2,:), '+') ; axis tight ; 
    pause ;
end

% state and switch vectors
xs = [xc; xf(:)] ;

%
%
%
function x = camera_matrix_to_pose(P)
x = [P(1:3,4); R2w(P(1:3,1:3))'] ;
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
function im = get_bundle_images(options, idx)
ncams = length(idx);
im = cell(ncams,1); % images cell
%warning('left images were cropped');
%fprintf('images: ');
for i = idx
    %fprintf(['#',num2str(i), ', ']);
    if mod(i,2);
        I = imread(strcat(options.img_folder{1},options.cam_left.image{(i+1)/2}));
        if size(I,3) == 1 % bayer decoding?
            I = demosaic(I, 'grbg');
        end
        %I = I(5:end-4,5:end-4,:);
    else
        I = imread(strcat(options.img_folder{2},options.cam_right.image{(i+0)/2}));
        if size(I,3) == 1 % bayer decoding?
            I = demosaic(I, 'grbg');
        end
    end
    if length(size(I)) == 3
        %I = I(:,:,2);
        I = rgb2gray(I);
    end
    im{i-idx(1)+1} = I;
end
%
function r = isdegenerate(x1,x2)
r = iscolinear(x1(1,:), x1(2,:), x2(1,:)) || ...
    iscolinear(x1(1,:), x1(2,:), x2(2,:)) || ...
    iscolinear(x1(1,:), x2(1,:), x2(2,:)) || ...
    iscolinear(x1(2,:), x2(1,:), x2(2,:));
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