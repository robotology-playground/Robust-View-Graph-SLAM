%function demo_graph_optimiser
% A demo of robust pose-graph constraints optimiser
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

%clear; 
clc; close;
%dbstop if error;
cd('/home/tabuhashim/Dropbox/code/matlab/icub');
run('./set_folders');

%folder = '~/Documents/data/blue/run_20160213';
%folder = '~/Documents/data/blue/run_20160224';
%folder = '~/Documents/data/blue/run_20160301';
%folder = '~/Documents/data/blue/run_20160319'; % FUSION16
%folder = '~/Documents/data/blue/run_20160327';
folder = '~/Documents/data/WalkingDatasets_1/results/M:20'; % walking
load(strcat(folder,'/options'));

options.save = '/home/tabuhashim/Documents/data/WalkingDatasets_1/results/M:20';

% get globals from kinematics
options.root = 0; % robot root, P{1} = A0;
load(strcat(options.save,'/encoders'));
try load(strcat(options.save, '/floatingbase'));
    [Pkin, A0, H0, A_left, A_right] = cameras_from_kinematics(eyes, neck, waist, floatingbase);
catch
    [Pkin, A0, H0, A_left, A_right] = cameras_from_kinematics(eyes, neck, waist);
end

% get constraints
load(strcat(options.save, '/constraints'));
SELECT_SOME_NODES = [0,0]; % IDs of first and last nodes to consider
if ~any(SELECT_SOME_NODES)==0 && length(SELECT_SOME_NODES)==2
    edges = vertcat(C.edge);
    nodes = unique(sort(edges(:)));
    f1 = max(SELECT_SOME_NODES(1), min(nodes));
    f2 = min(SELECT_SOME_NODES(2), max(nodes));
    C = get_constraints(C,f1,f2);
    scan = scan(f1:f2);
    kpts = kpts(f1:f2);
    Pkin = Pkin(f1:f2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get global poses : kinematics
M = size(Pkin, 2);
pose1 = zeros(6, M);
for i = 1:M
    pose1(:,i) = [Pkin{i}(1:3,4); R2w(Pkin{i}(1:3,1:3))'];
end
pose1 = transform_to_relative_w(pose1, pose1(:,1)); % move to 1st cam ref frame

close all;
figure(1);
subplot(3,1,1); plot(pose1(1,1:2:M),'b'); hold on; title('t - left');
subplot(3,1,2); plot(pose1(2,1:2:M),'b'); hold on;
subplot(3,1,3); plot(pose1(3,1:2:M),'b'); hold on; drawnow;

figure(2);
subplot(3,1,1); plot(pose1(4,1:2:M)*180/pi,'b'); hold on; title('a - left');
subplot(3,1,2); plot(pose1(5,1:2:M)*180/pi,'b'); hold on;
subplot(3,1,3); plot(pose1(6,1:2:M)*180/pi,'b'); hold on; drawnow;

figure(3);
plot3(pose1(1,1:2:M),pose1(2,1:2:M),pose1(3,1:2:M),'b'); hold on;
plot3(pose1(1,2:2:M),pose1(2,2:2:M),pose1(3,2:2:M),'b--');
if exist('floatingbase','var')
    plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:),'k');
end
drawnow; hold on; axis equal
% project cameras to floating-base
if 0
    pose1_f = pose1;
    for i=1:2:M % left camera
        H=(Pkin{i}/A_left{(i+1)/2})/H0;
        pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
    end
    for i=2:2:M % right camera
        H=(Pkin{i}/A_right{i/2})/H0;
        pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
    end
    figure(3);
    plot3(pose1_f(1,1:2:M),pose1_f(2,1:2:M),pose1_f(3,1:2:M),'r');hold on;
    plot3(pose1_f(1,2:2:M),pose1_f(2,2:2:M),pose1_f(3,2:2:M),'r--');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g2o by Rainer Kuemmerle; Giorgio Grisetti; Hauke Strasdat; Kurt Konolige
% and Wolfram Burgard;
%
% (insert code here)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotation averaging by Avishek Chatterjee and Venu Madhav Govindu
I = zeros(2,length(C));
RR = zeros(3,3,length(C));
for i = 1:length(C)
    I(:,i) = C(i).edge';
    RR(:,:,i) = w2R(C(i).z(4:6))';
end
R = AverageSO3Graph(RR, I, 'SIGMA', 5, 'MaxIterations', [1000, 1000]);
pose2 = pose1;
for i = 1:size(R,3)
    pose2(4:6,i) = R2w(R(:,:,i)')';
end
%pose2 = transform_to_global_w(pose2, pose1(:,1)); % move to kin ref frame

figure(2);
subplot(3,1,1); plot(pose2(4,1:2:end)*180/pi,'g'); hold on;
subplot(3,1,2); plot(pose2(5,1:2:end)*180/pi,'g'); hold on;
subplot(3,1,3); plot(pose2(6,1:2:end)*180/pi,'g'); hold on; drawnow;
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1DSfM by Kyle Wilson and	Noah Snavely
%
% (insert code here)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process constraints (resolve redundancies, and add weights)
C = assign_constraint_weight(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % get global poses : sequential
% pose = pose_generate_sequential(C2);
% figure;
% pose0 = pose;
% plot3(pose0(1,1:2:end),pose0(2,1:2:end),pose0(3,1:2:end),'b-*');
% hold on
% plot3(pose0(1,2:2:end),pose0(2,2:2:end),pose0(3,2:2:end),'r-*');
% axis equal; grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get global poses : maximum spanning tree
pose3 = pose_generate_spanning_tree(C);
%pose3 = transform_to_global_w(pose3, pose1(:,1)); % move to kin ref frame

figure(1);
subplot(3,1,1); plot(pose3(1,1:2:end),'m'); hold on;
subplot(3,1,2); plot(pose3(2,1:2:end),'m'); hold on;
subplot(3,1,3); plot(pose3(3,1:2:end),'m'); hold on; drawnow;
figure(2);
subplot(3,1,1); plot(pose3(4,1:2:end)*180/pi,'m'); hold on;
subplot(3,1,2); plot(pose3(5,1:2:end)*180/pi,'m'); hold on;
subplot(3,1,3); plot(pose3(6,1:2:end)*180/pi,'m'); hold on; drawnow;
figure(3);
plot3(pose3(1,1:2:end), pose3(3,1:2:end), -pose3(2,1:2:end), 'm');
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get global poses : optimiser
[~, ~, sw, x, C, Ct] = run_pose_graph_estimation(C);
pose4 = reshape(x, 6, length(x)/6);
%pose4 = transform_to_global_w(pose4, pose1(:,1)); % move to kin ref frame
%save(strcat(options.save,'/pose'),'sw', 'x', 'C', 'Ct');

figure(1);
subplot(3,1,1); plot(pose4(1,1:2:end),'r'); hold on;
subplot(3,1,2); plot(pose4(2,1:2:end),'r'); hold on;
subplot(3,1,3); plot(pose4(3,1:2:end),'r'); hold on; drawnow;

figure(2);
subplot(3,1,1); plot(pose4(4,1:2:end)*180/pi,'r'); hold on;
subplot(3,1,2); plot(pose4(5,1:2:end)*180/pi,'r'); hold on;
subplot(3,1,3); plot(pose4(6,1:2:end)*180/pi,'r'); hold on; drawnow;

figure(3);
plot3(pose4(1,1:2:end), pose4(3,1:2:end), -pose4(2,1:2:end), 'r');

figure(4); plot_graph_constraints(pose4, C, sw, .03, .025);
xlabel('Eastern (m)'); ylabel('Depth (m)'); box on;
axis equal; view(2); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare parameters
figure(5);
% for i = 1:3
%     subplot(3,2,2*(i-1)+1); hold on;   % left camera translation
%     plot(pose1(i,1:2:end),'g');
%     plot(pose2(i,1:2:end),'r');
%     plot(pose3(i,1:2:end),'b');
%     plot(pose4(i,1:2:end),'k');
% end
for i = 4:6
    subplot(3,1,i-3);
    %subplot(3,2,2*(i-3)-1);   % left camera
    hold on;
    %plot((pose1(i,1:2:end))*180/pi,'k'); % kinematics
    plot((pose2(i,1:2:end)-pose1(i,1:2:end))*180/pi,'g'); % averaging
    plot((pose3(i,1:2:end)-pose1(i,1:2:end))*180/pi,'m'); % mst
    plot((pose4(i,1:2:end)-pose1(i,1:2:end))*180/pi,'r'); % RNLS
    %     subplot(3,2,2*(i-3));  % right camera
    %     hold on;
    %     %plot((pose1(i,2:2:end))*180/pi,'k');
    %     plot((pose2(i,2:2:end)-pose1(i,1:2:end))*180/pi,'g');
    %     plot((pose3(i,2:2:end)-pose1(i,1:2:end))*180/pi,'m');
    %     plot((pose4(i,2:2:end)-pose1(i,1:2:end))*180/pi,'r');
end
title('errors wrt kinematics');

% stereo angles
stereo_pose_k=zeros(6,size(pose1,2)/2);
stereo_pose_a=zeros(6,size(pose1,2)/2);
stereo_pose_m=zeros(6,size(pose1,2)/2);
stereo_pose_n=zeros(6,size(pose1,2)/2);
for i=1:size(pose1,2)/2
    stereo_pose_k(:,i) = transform_to_relative_w(pose1(:,2*i), pose1(:,2*i-1)); % kinematics
    stereo_pose_a(:,i) = transform_to_relative_w(pose2(:,2*i), pose2(:,2*i-1)); % averaging
    stereo_pose_m(:,i) = transform_to_relative_w(pose3(:,2*i), pose3(:,2*i-1)); % MST
    stereo_pose_n(:,i) = transform_to_relative_w(pose4(:,2*i), pose4(:,2*i-1)); % RNLS
end
figure(6);
subplot(2,2,1); 
plot([eyes; neck]'*180/pi,'linewidth', 2);grid on;
legend('tilt','pan','vergence','head roll','head pitch','head yaw');
subplot(2,2,2); 
plot(stereo_pose_k(4,:)*180/pi,'b'); hold on;
plot(stereo_pose_a(4,:)*180/pi,'g');
plot(stereo_pose_m(4,:)*180/pi,'m');
plot(stereo_pose_n(4,:)*180/pi,'r');
subplot(2,2,3); 
plot(stereo_pose_k(5,:)*180/pi,'b'); hold on;
plot(stereo_pose_a(5,:)*180/pi,'g');
plot(stereo_pose_m(5,:)*180/pi,'m');
plot(stereo_pose_n(5,:)*180/pi,'r');
subplot(2,2,4); 
plot(stereo_pose_k(6,:)*180/pi,'b'); hold on;
plot(stereo_pose_a(6,:)*180/pi,'g');
plot(stereo_pose_m(6,:)*180/pi,'m');
plot(stereo_pose_n(6,:)*180/pi,'r');

figure(7);
subplot(3,1,1); 
plot(stereo_pose_k(1,:),'b'); hold on;
plot(stereo_pose_m(1,:),'m');
plot(stereo_pose_n(1,:),'r');
subplot(3,1,2); 
plot(stereo_pose_k(2,:),'b'); hold on;
plot(stereo_pose_m(2,:),'m');
plot(stereo_pose_n(2,:),'r');
subplot(3,1,3); 
plot(stereo_pose_k(3,:),'b'); hold on;
plot(stereo_pose_m(3,:),'m');
plot(stereo_pose_n(3,:),'r');

%[options,encoders]=set_images(options);
%figure(8);
%plot( options.cam_right.time , 'linewidth', 2);

%%
figure(9); clf;
if exist('floatingbase','var')
    plot3(-floatingbase(2,:),floatingbase(1,:),floatingbase(3,:),'k');hold on;
end
% % project (MST) cameras to floating-base
% ref_pose = [Pkin{1}(1:3,4); R2w(Pkin{1}(1:3,1:3))'];
% pose1_f = transform_to_global_w(pose3, ref_pose); % move to kin ref frame
% M = size(pose1_f,2);
% for i=1:2:M % left camera
%     A=[w2R(pose1_f(4:6,i)), pose1_f(1:3,i)];
%     H=(A/A_left{(i+1)/2})/H0;
%     pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
% end
% for i=2:2:M % right camera
%     A=[w2R(pose1_f(4:6,i)), pose1_f(1:3,i)];
%     H=(A/A_right{i/2})/H0;
%     pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
% end
% plot3(-pose1_f(2,1:2:M),pose1_f(1,1:2:M),pose1_f(3,1:2:M),'b');
% plot3(-pose1_f(2,2:2:M),pose1_f(1,2:2:M),pose1_f(3,2:2:M),'b--');
% project (graph) cameras to floating-base
ref_pose = [Pkin{1}(1:3,4); R2w(Pkin{1}(1:3,1:3))'];
pose1_f = transform_to_global_w(pose4, ref_pose); % move to kin ref frame
M = size(pose1_f,2);
for i=1:2:M % left camera
    A=[w2R(pose1_f(4:6,i)), pose1_f(1:3,i)];
    H=(A/A_left{(i+1)/2})/H0;
    pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
end
for i=2:2:M % right camera
    A=[w2R(pose1_f(4:6,i)), pose1_f(1:3,i)];
    H=(A/A_right{i/2})/H0;
    pose1_f(:,i)=[H(1:3,4);R2w(H(1:3,1:3))'];
end
plot3(-pose1_f(2,1:2:M),pose1_f(1,1:2:M),pose1_f(3,1:2:M),'r'); hold on;
plot3(-pose1_f(2,2:2:M),pose1_f(1,2:2:M),pose1_f(3,2:2:M),'r--');

xlabel('Eastern (m)'); ylabel('Depth (m)'); box on;
axis equal; view(2); grid on; drawnow;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get global poses w.r.t. icub root
pose5 = pose4;
P = cell(1, size(pose5, 2));
for i = 1:size(pose5, 2)
    %M = [w2R(pose(4:6,i)), pose(1:3,i); 0 0 0 1];
    M = [A0(1:3,1:3) A0(1:3,4)/1000; 0 0 0 1]*[w2R(pose5(4:6,i)), pose5(1:3,i); 0 0 0 1];
    P{i} = M(1:3,1:4);
    %P{i}=[w2R(pose4(4:6,i)), pose4(1:3,i)];
    pose5(1:3, i) = P{i}(1:3, 4);
    pose5(4:6, i) = R2w(P{i}(1:3, 1:3))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove outlier scans
use = zeros(1,size(pose5, 2));
for i = 1:length(C);
    %if C(i).edge(2)-C(i).edge(1) == 1 && mod(C(i).edge(1),2) == 1
    if abs(C(i).edge(2)-C(i).edge(1)) == 1
        if sw(i) == 1;
            use(C(i).edge(1)) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODIFY_MIN_VIEWS = 0;
if MODIFY_MIN_VIEWS % remove more scan outliers?
    nviews = min(MODIFY_MIN_VIEWS,options.ncams);
    [scan, kpts] = remove_more_scan_outliers(options, pose, C, nviews, scan, kpts, ncam);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot map and trajectory
figure(7);
%check_depth_map(scan(k), kpts(k), options, P(k), 1200, 1);
check_depth_map(scan, kpts, options, Pkin, 2, 0, use);
figure(8); %hold on;
check_depth_map(scan, kpts, options, P, 2, 0, use);
%plot3(-pose5(3,1:2:end),pose5(2,1:2:end),pose5(1,1:2:end),'g','linewidth', 2);
plot3(-pose5(2,1:2:end),pose5(1,1:2:end),pose5(3,1:2:end),'g','linewidth', 2);
hold on
%plot3(-pose5(3,2:2:end),pose5(2,2:2:end),pose5(1,2:2:end),'m','linewidth', 2);
plot3(-pose5(2,1:2:end),pose5(1,1:2:end),pose5(3,1:2:end),'g','linewidth', 2);
%view(-31.6000,32.4000);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dense mapping using graph estimate

pose5 = pose4;
%pose5(4:6,:) = pose2(4:6,:);

% P = cell(1, size(pose5, 2));
% for i = 1:size(pose5, 2)
%     %M = [w2R(pose(4:6,i)), pose(1:3,i); 0 0 0 1];
%     M = [w2R(pose5(4:6,i)), pose5(1:3,i); 0 0 0 1];
%     P{i} = M(1:3,1:4);
%     %P{i}=[w2R(pose4(4:6,i)), pose4(1:3,i)];
%     pose5(1:3, i) = P{i}(1:3, 4);
%     pose5(4:6, i) = R2w(P{i}(1:3, 1:3))';
% end
%pose5 = transform_to_global_w(pose4, pose1(:,1)); % move to kin ref frame
P = cell(1, size(pose5, 2));
for i = 1:size(pose5, 2)
    M = [w2R(pose5(4:6,i)), pose5(1:3,i); 0 0 0 1];
    P{i} = M(1:3,1:4);
end

for i=459;%:2:size(pose5,2)
    i
    depth_sat = 6;
    im1=imread(options.cam_left.image{(i+1)/2});
    im2=imread(options.cam_right.image{(i+1)/2});
    % convert images to gray scale
    if size(im1,3) > 1; im1 = rgb2gray(im1); end;
    if size(im2,3) > 1; im2 = rgb2gray(im2); end;
    p=transform_to_relative_w(pose5(:,i+1),pose5(:,i));
    t=p(1:3);
    R=w2R(p(4:6));
    Sv=cv.stereoRectify(options.K1,options.kc1,options.K2,options.kc2,...
        size(im1'),R', t, 'Alpha',0,'ZeroDisparity',1);
    [mapx1v, mapy1v]=cv.initUndistortRectifyMap(...
        options.K1,options.kc1,Sv.P1(1:3,1:3),size(im1'),'R',Sv.R1);
    [mapx2v, mapy2v]=cv.initUndistortRectifyMap(...
        options.K2,options.kc2,Sv.P2(1:3,1:3),size(im2'),'R',Sv.R2);
    dst1v=cv.remap(im1,mapx1v,mapy1v);
    dst2v=cv.remap(im2,mapx2v,mapy2v);
    Q = Sv.Q;
    my_params();
    [D1, D2] = elasMex(dst1v', dst2v', param);
    
    figure(8); clf; imshow(im1);
    %subplot(1,3,1); %imagesc(D1');
    
    % project disparity to 3d
    im3d = cv.reprojectImageTo3D(D1, Q);
    X = -im3d(:,:,2);
    Y =  im3d(:,:,1);
    Z = -im3d(:,:,3);
    
    % saturate depth and remove bad disparity.
    depth = sqrt(X.^2+Y.^2+Z.^2);
    index = (Z>0 & depth<depth_sat);
    points_3d=[];
    X = X(index); points_3d(:,1) = X(:);
    Y = Y(index); points_3d(:,3) = Y(:);
    Z = Z(index); points_3d(:,2) = Z(:);
    if isempty(points_3d); continue; end
    
    % if you use pose5. everything is going to be in camera
    % coordinates.
    
    % if you use pose1_f. everything is going to be in floating base
    % coordinates.
    
    ref_pose = [Pkin{1}(1:3,4); R2w(Pkin{1}(1:3,1:3))'];
    pose5_g=transform_to_global_w(pose5,ref_pose);
    
    points_3d=transform_to_global_w(points_3d,pose1_f(:,i));
    figure(9); clf;
    plot3(points_3d(:,1),points_3d(:,2),points_3d(:,3),'.','markersize',.5);
    hold on
    
    %plot3(-pose5_g(2,1:2:end),pose5_g(1,1:2:end),pose5_g(3,1:2:end),'g','linewidth', 2);
    %plot3(-pose5_g(2,2:2:end),pose5_g(1,2:2:end),pose5_g(3,2:2:end),'m','linewidth', 2);
    
    plot3(-pose1_f(2,1:2:end),pose1_f(1,1:2:end),pose1_f(3,1:2:end),'g','linewidth', 2);
    plot3(-pose1_f(2,2:2:end),pose1_f(1,2:2:end),pose1_f(3,2:2:end),'m','linewidth', 2);
    plot3(-floatingbase(2,:),floatingbase(1,:),floatingbase(3,:),'k','linewidth', 2);
    axis equal; view(-36.6,38.8); drawnow;
    
    %figure(9); clf;
    %use=zeros(1,size(P,2)); use(i)=1;
    %check_depth_map(scan, kpts, options, P, 5, 0, use, 1);
    %hold on
    %plot3(pose5(1,1:2:end),pose5(3,1:2:end),-pose5(2,1:2:end),'g','linewidth', 2);
    %plot3(pose5(1,2:2:end),pose5(3,2:2:end),-pose5(2,2:2:end),'m','linewidth', 2);
    %axis equal;
    
    %pause;
end