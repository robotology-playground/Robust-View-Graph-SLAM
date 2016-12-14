clear; clc;
%run('/home/tabuhashim/Dropbox/code/matlab/icub/set_folders');

folder=['~/Documents/data/WalkingDatasets_1/results','/run_20160804'];
load(strcat(folder,'/options'));
load(strcat(folder,'/bundle_1'));
load(strcat(folder,'/C_1'));

% Data
ncams = size(p,2);
npnts = sum(p(1).s);
ncnst = length(C);
fileID = fopen('/home/tabuhashim/problem-12-pre.txt','w');
fprintf(fileID,'%d %d %d\n',[ncams,npnts,ncnst]);
for i = 1 : ncnst
    data(i,:) = [C(i).cam-1 C(i).kpt-1 C(i).z'];
end
% measurements
for i = 1 : size(data,1)
    fprintf(fileID,'%d %d     %f %f\n',data(i,:));
end
% cams
for i = 1 : ncams
    xc = [xs((i-1)*6+(4:6)); xs((i-1)*6+(1:3))];
    %xc = transform_to_relative_w(zeros(6,1),xc);
    fprintf(fileID,'%f\n%f\n%f\n%f\n%f\n%f\n',xc);
    fprintf(fileID,'%f\n%f\n%f\n',[1 0 0]);
end
% 3D points
for i = 1 : npnts
    if p(1).s(i)==0
        continue
    end
    q = [p(1).x(i);p(1).y(i)];
    r = 1/xs((6*ncams)+i);
    rim = sqrt(q(1)*q(1) + q(2)*q(2) + 1);
    d = r./rim;
    Q(1) = d*q(1);
    Q(2) = d*q(2);
    Q(3) = d;
    fprintf(fileID,'%f\n%f\n%f\n',Q);
end
fclose(fileID);

% Ceres
exect = '/home/tabuhashim/Dev/ceres-bin/bin/bundle_adjuster';
input = '--input=/home/tabuhashim/problem-12-pre.txt';
output = '-final_ply /home/tabuhashim/problem-12-post.txt';
params = ' ';%'-translation_sigma .000001 -rotation_sigma .01 -point_sigma .01';
options = '-num_iterations 1000 -robustify';
tic
status = system([exect ' ' input ' ' output ' ' params ' ' options]);
toc

% a = load('/home/tabuhashim/problem-12-post.txt');
% pose1 = xs(1:6);
% pose2 = xs(7:12);
% init = transform_to_relative_w(pose2, pose1)
% pose1 = a(1,1:6)';
% pose2 = a(2,1:6)';
% est = transform_to_relative_w(pose2, pose1)

load(strcat(folder,'/constraints'));
pose_batch=zeros(7,length(C)+1);
pose_ceres=zeros(7,length(C)+1);
a = load('/home/tabuhashim/Output.txt');
pose = a(:,[4:6 1:3])';
for i=1:length(C)
    cam = mod(C(i).edge(2),2) + 1; % 0 is for no camera, 1 for even, 2 for odd
    pose_batch(:,C(i).edge(2)) = [C(i).z; cam];
    pose_ceres(:,C(i).edge(2)) = [transform_to_relative_w(pose(:,i+1), pose(:,1)); cam];
end

cam = pose_batch(7,:);

% figure;
% subplot(3,1,1); plot(pose_batch(4,cam==2)*180/pi,'r');
% hold on; plot(pose_ceres(4,cam==2)*180/pi,'r--');
% subplot(3,1,2); plot(pose_batch(5,cam==2)*180/pi,'r');
% hold on; plot(pose_ceres(5,cam==2)*180/pi,'r--');
% subplot(3,1,3); plot(pose_batch(6,cam==2)*180/pi,'r');
% hold on; plot(pose_ceres(6,cam==2)*180/pi,'r--');

figure;
subplot(3,1,1); plot(pose_batch(4,cam==1)*180/pi,'r');
hold on; plot(pose_ceres(4,cam==1)*180/pi,'b');
subplot(3,1,2); plot(pose_batch(5,cam==1)*180/pi,'r');
hold on; plot(pose_ceres(5,cam==1)*180/pi,'b');
subplot(3,1,3); plot(pose_batch(6,cam==1)*180/pi,'r');
hold on; plot(pose_ceres(6,cam==1)*180/pi,'b');

% figure;
% subplot(3,1,1); plot(pose_batch(1,cam==2),'r');
% hold on; plot(pose_ceres(1,cam==2),'b');
% subplot(3,1,2); plot(pose_batch(2,cam==2),'r');
% hold on; plot(pose_ceres(2,cam==2),'b');
% subplot(3,1,3); plot(pose_batch(3,cam==2),'r');
% hold on; plot(pose_ceres(3,cam==2),'b');

figure;
subplot(3,1,1); plot(pose_batch(1,cam==1),'r');
hold on; plot(pose_ceres(1,cam==1),'b');
subplot(3,1,2); plot(pose_batch(2,cam==1),'r');
hold on; plot(pose_ceres(2,cam==1),'b');
subplot(3,1,3); plot(pose_batch(3,cam==1),'r');
hold on; plot(pose_ceres(3,cam==1),'b');