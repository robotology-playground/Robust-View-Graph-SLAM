function [C, scan] = batch_motion_and_map_inverse_depth_v2(options, Pkin)
%[C, scan] = batch_motion_and_map_inverse_depth_v2(options, P)
% Robust epipolar constraints optimiser
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot
%

% Configurations
switch_config;

% Adjust ncams for last few image frames
nimages = 2*size(options.cam_left.image, 2); % left + right = 2*left
ncams = min(options.ncams, nimages);
%nkeys = min(options.nkeys, ncams);
scan = cell(1, ncams);
kpts = cell(1, ncams);
ncam = cell(1, ncams); % number of cams seeing each point in each bundle
npts = cell(1, ncams); % number of points tracked in each cam in each bundle
C = struct('edge', [], 'z', [], 'c', []); % Graph constraints, empty structure fields

% Check for previous computations ;
% for instance, if code was not completed, then continue ;
count = 0 ; % counter for number of already computed constraints (used later)
if exist(strcat(options.save,'/constraints.mat'), 'file') == 2
    [C, scan, kpts, ncam, npts] = load_constraints(options);
end

k = 1; % the first input image

while k <  nimages;%nkeys+1
    
    % Find (new) ncams, nkeys and image indices
    ncams = min(nimages-k+1, options.ncams);
    nkeys = min(options.nkeys, ncams);
    if k <= size(scan, 2) % check if this was done before for this camera
        if ~isempty(scan{k});
            %k = k + round(ncams/nkeys) - mod(round(ncams/nkeys),2);
            k = k + ceil(ncams/nkeys);
            continue
        end
    end
    
    % Determine a geometrically stable ncams
%     if ~isempty(Pkin)
%         A = [Pkin{k}; 0 0 0 1]; i=0; t=[];
%         while (k+i)<nimages
%             M=A\[Pkin{k+i}; 0 0 0 1];
%             t(i+1)=abs(M(1,4));
%             if sum(t>.05)>(0.75*ncams) % 75% of poses
%                 ncams=i;
%                 break
%             end
%             i=i+1;
%         end
%     end
    
    % Get conrners and tracks relative to the reference image
    % Also, update ncams accordingly and then compute npnts (num. of points)
    p = get_aligned_point_matches_v2(options, k, ncams);
    ncams = length(p);% number of images (in case some were removed)
    npnts = sum(p(1).s);% number of points tracked, the same for all images
    fprintf('Data contains %d image and %d tracks.\n',ncams,npnts);
    
    % Get the kinematics
    if ~isempty(Pkin)
        idx = zeros(1,ncams);
        for i = 1:ncams
            idx(i) = p(i).c;
        end
        M = Pkin(idx);
    else
        M = [];
    end
    if 0; 
        figure;
        plot_5(Pkin, idx); 
        view(2); pause
    end

    % First, the state vector is initialised
    % using kinematics and vision. It is considered as the only main input
    % to the PwgOptimiser class;
    [xs, p] = initialise_linearisation_point_v2(options, p, M);
    ncams = length(p);% number of images (in case some were removed)
    npnts = sum(p(1).s);% number of points tracked, the same for all images
    fprintf('Bundle contains %d image and %d tracks.\n',ncams,npnts);
    save(strcat(options.save,['/bundle_',num2str(k)]),'-v7.3','p','xs'); % exporting data to Ceres
    
    if 0;
        figure;
        plot_1(xs, p, ncams, M{1});
        pause;
    end
    
    return
    
    % Constraints are pushed into the PwgOptimiser class. These constraints were
    % initialised only from vision tracks;
    if size(p,2)<2 || sum(p(1).s)<10
        disp('Not enough data .......');
        continue
    end
    [CC, CCt, sw] = initialise_constraints(p, options);
    save(strcat(options.save,['/C_',num2str(p(1).c)]),'-v7.3','CC','CCt','sw'); % exporting data to Ceres
    
    % Optimise image constraints
    %options.fid = fopen([strcat(options.save) '/bundle_',num2str(k),'.txt'],'w');  % Open file
    [x, ~, CC, sw] = PwgOptimiser(CC, CCt, xs, sw, ncams, options);
    %[x, ~, CC, sw] = PwgOptimiser(xs, p, options);
    %fclose(options.fid); % Close file
    
    if 0;
        figure;
        plot_1(x, p, ncams, M{1}); pause;
    end
    
    % Generate scans, visibility, and remove unwanted points
    %nviews = min(ceil(min(options.ncams/2, ncams)), max(ncams/2, 2)); % at least 3 frames
    nviews = options.nview;
    vis = zeros(1, npnts);
    cam = zeros(1, ncams);
    for i = 1 : length(CC)
        vis(CC(i).kpt) = vis(CC(i).kpt) + sw(i); % number of cams seeing a point
        cam(CC(i).cam) = cam(CC(i).cam) + sw(i); % number of points tracked in each cam
    end
    %nviews = round(mean(vis))-1;
    %index{k} = find(vis>nviews);
    index = find(vis>=nviews);
    kpts{k} = [p(1).x(index); p(1).y(index)];%p{1}(:, index);
    scan{k} = x(ncams*6 + index)';
    ncam{k} = vis(index); % number of cams seeing a point in this bundle
    npts{k} = cam; % number of points tracked in each cam in this bundle
    if 0;
        clf; plot(x((ncams*6+1):end)); hold on;
        plot(index,x(ncams*6+index),'rO'); pause;
    end;
    
    % Generate graph constraints and scans
    for c = 2 : ncams
        count = count + 1;
        C(count).edge = [p(1).c, p(c).c];
        %C(count).t = x(c*6+(1:3));
        %C(count).a = x(c*6+(4:6));
        C(count).z = x((c-1)*6+(1:6));
        C(count).c = cam(c); % number of points tracked in each cam in
        % this bundle this can be used as an edge weight in the pose graph
        %C(count).Y=Y(c*6+(1:6),c*6+(1:6));
        %C(count).P=P(c*6+(1:6),c*6+(1:6));
    end
    
    % Display relative poses
    %for j = 1 : ncams
    %    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d\n', ...
    %        x((j-1)*6+(1:3))',x((j-1)*6+(4:6))'*180/pi,cam(j));
    %end
    
    % Save
    save(strcat(options.save,'/constraints'),'-v7.3', ...
        'C','kpts','scan','ncam','npts','USE_VISION');
    
    % Next reference frame
    %if (round(ncams/nkeys)-mod(round(ncams/nkeys),2)) == 0, break, end
    %k = k + round(ncams/nkeys) - mod(round(ncams/nkeys),2);
    k = k + ceil(options.ncams/options.nkeys); % next key-frame
    
    % Show map results
    if 0 ;%mod(k,2);
        clf; plot_1(x, p, ncams);
    end
    
    % Show map results
    if 0
        clf; plot_2(scan, ncams, npts, C, Pkin, kpts, options);
        pause;
    end
    
    % Show trajectory results
    if 0;
        clf; plot_3(x, Pkin, ncams);
    end
    
    % Display the outliers
    if 0
        clf; plot_4(k, options, CC);
    end
    
end

%
%
function plot_1(x, p, ncams, H)
if nargin<3
    H=eye(3,4);
end
ref_pose = camera_matrix_to_pose(H);        
r = 1./x((ncams*6+1):length(x))';
m1 = [p(1).x; p(1).y];
idx = r > 0; % remove negative depth
X = get_scan_from_range(m1(:,idx), r(idx));
X = transform_to_global_w(X, ref_pose);
plot3(X(2,:), X(1,:), X(3,:), '+') ; hold on;
cam = ones(1, ncams);
for i = 1 : ncams
    cam(i) = mod(p(i).c, 2);
end
pose = transform_to_global_w(reshape(x(1:6*ncams),6,ncams), ref_pose);
plot3(pose(2,cam==1), pose(1,cam==1), pose(3,cam==1), 'r', 'linewidth', 2);
plot3(pose(2,cam==0), pose(1,cam==0), pose(3,cam==0), 'g', 'linewidth', 2);
axis equal ; view(2); drawnow;
%
%
function plot_3(x,Pkin,ncams)
plot3(x(1:12:6*ncams),x(2:12:6*ncams),x(3:12:6*ncams),'b'); % left camera
hold on;plot3(x(7:12:6*ncams),x(8:12:6*ncams),x(9:12:6*ncams),'r'); % right camera
xk=camera_matrix_to_pose(Pkin);
plot3(xk(1:12:6*ncams),xk(2:12:6*ncams),xk(3:12:6*ncams),'b'); % left camera kinematics
plot3(xk(7:12:6*ncams),xk(8:12:6*ncams),xk(9:12:6*ncams),'r'); % right camera kinematics
pause;
%
%
function plot_4(k,options,CC)
im = get_bundle_images(options, k);
for j = 1 : ncams
    figure(j);
    imshow(im{j}); hold on;
end
for j = 1 : length(CC)
    kpt = CC(j).kpt;
    cam = CC(j).cam;
    [K, kc] = get_intrinsics(options, cam+k-1);
    [x, y] = UndistPointFor(p{cam}(1,kpt), p{cam}(2,kpt),...
        kc(1), kc(2), kc(3), kc(4), kc(5));
    x = K(1,1)*x + K(1,3);
    y = K(2,2)*y + K(2,3);
    figure(cam);
    if sw(j)==1
        plot(x, y, 'go');
    else
        plot(x, y, 'ro');
    end
end
pause;
%
%
function plot_5(Pkin,k)
%pose=zeros(6,size(Pkin,2));
%for i = 1 : size(Pkin,2)
%    pose(1:3,i)=Pkin{i}(1:3,4);
%    pose(4:6,i)=R2w(Pkin{i}(1:3,1:3))';
%end
pose = camera_matrix_to_pose(Pkin);
plot3(pose(2,1:2:end),pose(1,1:2:end),pose(3,1:2:end),'b+'); % all cams
hold on;
plot3(pose(2,2:2:end),pose(1,2:2:end),pose(3,2:2:end),'g+'); % all cams
plot3(pose(2,k),pose(1,k),pose(3,k),'rO'); % processed cameras
axis equal; axis tight; grid on; box on;
xlabel('Eastern (m)'); ylabel('Forward (m)');
drawnow;
%
%
function p=get_scan_from_range(p,r)
rim=sqrt(p(1,:).*p(1,:)+p(2,:).*p(2,:)+1);
d=r./rim;
p(1,:)=d.*p(1,:);
p(2,:)=d.*p(2,:);
p(3,:)=d;
%
%
function x=camera_matrix_to_pose(P)
if iscell(P);
    n=length(P);
    x=zeros(6,n);
    for i=1:n
        x(:,i)=[P{i}(1:3,4);R2w(P{i}(1:3,1:3))'];
    end
else
    x=[P(1:3,4);R2w(P(1:3,1:3))'];
end