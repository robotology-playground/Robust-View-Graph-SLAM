%function demo_constraints_optimiser
clear; clc; close all;
%dbstop if error;
cd('/home/tabuhashim/Dropbox/code/matlab/icub');
run('./set_folders');

NEW_OPTIONS = 1;
PROCESS_CONSTRAINTS = 1;
SELECT_SOME_NODES = 0; % IDs of first and last nodes to consider
MODIFY_MIN_VIEWS = 0;

if  NEW_OPTIONS
    robot = '/home/tabuhashim/Documents/data/blue';
    %options.folder = [robot,'/20150515_1'];
    %options.folder = [robot,'/20150901_3'];
    %options.folder = [robot,'/20151207_desk3'];
    options.folder = [robot,'/20151207_moving_2'];
    options.save = [robot,'/run_',datestr(now,'yyyymmdd')];
    options.first_image	= 101;%51;
    options.last_image = 2100;%550;
    options.steps = 5;
    options.freq = 15;
    options = set_params(options);
    [options, encoders] = set_images(options);
else
    %folder = '/home/tabuhashim/Documents/data/blue/run_20160213';
    %folder = '/home/tabuhashim/Documents/data/blue/run_20160224';
    %folder = '/home/tabuhashim/Documents/data/blue/run_20160301';
    %folder = '/home/tabuhashim/Documents/data/blue/run_20160319';
    folder = '/home/tabuhashim/Documents/data/blue/run_20160327';
    load(strcat(folder,'/options'));
end

% get globals from kinematics
options.root = 0; % robot root, P{1} = A0;
load(strcat(options.save,'/encoders'));
[Pkin, A0] = cameras_from_kinematics(eyes, neck, waist);

if PROCESS_CONSTRAINTS
    % 2-views implementation
    %[kpts, desc] = batch_features(options);
    %[C, G] = batch_matching_new(kpts, desc, options);
    %[C, P] = batch_constraints_from_kinematics_new(encoders, C, options);
    %[C, scan] = batch_Emat_v2_new(C, kpts, options, P);
    
    % M-views implementation
    [C, scan] = batch_motion_and_map_inverse_depth_v2(options, Pkin);
end
load(strcat(options.save, '/constraints'));

if ~any(SELECT_SOME_NODES)==0 && length(SELECT_SOME_NODES)==2 && ~PROCESS_CONSTRAINTS
    edges = vertcat(C.edge);
    nodes = unique(sort(edges(:)));
    f1 = max(SELECT_SOME_NODES(1), min(nodes));
    f2 = min(SELECT_SOME_NODES(2), max(nodes));
    C = get_constraints(C, f1, f2);
    scan = scan(f1:f2);
    kpts = kpts(f1:f2);
    Pkin = Pkin(f1:f2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process constraints (resolve redundancies, and add weights)
C = assign_constraint_weight(C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get global poses : kinematics
pose1 = zeros(6, size(Pkin,2));
for i = 1:size(pose1, 2)
    M = [Pkin{1}; 0 0 0 1]\[Pkin{i}; 0 0 0 1];
    pose1(1:3,i) = M(1:3,4);
    pose1(4:6,i) = R2w(M(1:3,1:3))';
end

figure;
subplot(3,1,1); plot(pose1(4,1:2:end)*180/pi, 'b'); hold on;
subplot(3,1,2); plot(pose1(5,1:2:end)*180/pi, 'b'); hold on;
subplot(3,1,3); plot(pose1(6,1:2:end)*180/pi, 'b'); hold on;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g2o by Rainer Kuemmerle; Giorgio Grisetti; Hauke Strasdat; Kurt Konolige
% and Wolfram Burgard;
%
% (insert code here)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation averaging by Avishek Chatterjee and Venu Madhav Govindu
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
subplot(3,1,1); plot(pose2(4,1:2:end)*180/pi, 'g'); hold on;
subplot(3,1,2); plot(pose2(5,1:2:end)*180/pi, 'g'); hold on;
subplot(3,1,3); plot(pose2(6,1:2:end)*180/pi, 'g'); hold on;
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1DSfM by Kyle Wilson and	Noah Snavely
%
% (insert code here)
%


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
subplot(3,1,1); plot(pose3(4,1:2:end)*180/pi, 'k'); hold on;
subplot(3,1,2); plot(pose3(5,1:2:end)*180/pi, 'k'); hold on;
subplot(3,1,3); plot(pose3(6,1:2:end)*180/pi, 'k'); hold on;
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get global poses : optimiser
[y, Y, sw, x, C, Ct] = run_pose_graph_estimation(C);
pose4 = reshape(x, 6, length(x)/6);
save(strcat(options.save,'/pose'),'sw', 'x', 'C', 'Ct');
figure;
plot_graph_constraints(pose4, C, sw, .005, .01); view(2);
%xlabel('Eastern (m)'); ylabel('Depth (m)'); box on;
% figure;
% plot3(pose2(1,1:2:end),pose2(3,1:2:end),pose2(2,1:2:end),'g-x');
% hold on
% plot3(pose2(1,2:2:end),pose2(3,2:2:end),pose2(2,2:2:end),'m-x');
% axis equal; grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare parameters
figure;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get global poses w.r.t. icub root
pose5 = pose4;
P = cell(1, size(pose5, 2));
for i = 1:size(pose5, 2)
    %M = [w2R(pose(4:6,i)), pose(1:3,i); 0 0 0 1];
    M = [A0(1:3,1:3) A0(1:3,4)/1000; 0 0 0 1]*[w2R(pose5(4:6,i)), pose5(1:3,i); 0 0 0 1];
    P{i} = M(1:3,1:4);
    pose5(1:3, i) = P{i}(1:3, 4);
    pose5(4:6, i) = R2w(P{i}(1:3, 1:3))';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove outlier scans
use = zeros(1,size(pose5, 2));
for i = 1:length(C);
    %if C(i).edge(2)-C(i).edge(1) == 1 && mod(C(i).edge(1),2) == 1
    if abs(C(i).edge(2)-C(i).edge(1)) == 1
        if sw(i) == 1;
            use(C(i).edge(1)) = 1;
        end
    end
end
figure;
%check_depth_map(scan(k), kpts(k), options, P(k), 1200, 1);
check_depth_map(scan, kpts, options, Pkin, 1.5, 0, 0, use);
figure; %hold on;
check_depth_map(scan, kpts, options, P, 10, 0, 0, use);
plot3(-pose5(3,1:2:end),pose5(2,1:2:end),pose5(1,1:2:end),'g','linewidth', 2);
hold on
plot3(-pose5(3,2:2:end),pose5(2,2:2:end),pose5(1,2:2:end),'m','linewidth', 2);
%view(-31.6000,32.4000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MODIFY_MIN_VIEWS && ~PROCESS_CONSTRAINTS % remove more scan outliers?
    nviews = min(MODIFY_MIN_VIEWS,options.ncams);
    [scan, kpts] = remove_more_scan_outliers(options, nviews, C, scan, kpts, Pkin, ncam);
end
