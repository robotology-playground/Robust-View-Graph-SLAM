function plot_pose_parameters(pose,style,handle)

left_pose=pose(:,1:2:end);
right_pose=pose(:,2:2:end);
xlimit=size(pose,2)/2;

% translation left_pose
if nargin<4;figure(1);else figure(handle(1)); end;
subplot(3,1,1);plot(left_pose(1,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('x-dir (mm.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,2);plot(left_pose(2,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('y-dir (mm.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,3);plot(left_pose(3,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('z-dir (mm.)');
axis([0,xlimit,-inf,inf]);
xlabel('time steps');

% rotation left_pose
if nargin<4;figure(2);else figure(handle(2)); end;
subplot(3,1,1);plot(left_pose(4,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('tilt (deg.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,2);plot(left_pose(5,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('heading (deg.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,3);plot(left_pose(6,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('pan (deg.)');
axis([0,xlimit,-inf,inf]);
xlabel('time steps');

% translation right_pose
if nargin<4;figure(3);else figure(handle(3)); end;
subplot(3,1,1);plot(right_pose(1,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('x-dir (mm.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,2);plot(right_pose(2,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('y-dir (mm.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,3);plot(right_pose(3,:),style{1},'color',style{2});
grid on,box on,hold on,ylabel('z-dir (mm.)');
axis([0,xlimit,-inf,inf]);
xlabel('time steps');

% rotation right_pose
if nargin<4;figure(4);else figure(handle(4)); end;
subplot(3,1,1);plot(right_pose(4,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('tilt (deg.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,2);plot(right_pose(5,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('heading (deg.)');
axis([0,xlimit,-inf,inf]);
subplot(3,1,3);plot(right_pose(6,:)*180/pi,style{1},'color',style{2});
grid on,box on,hold on,ylabel('pan (deg.)');
axis([0,xlimit,-inf,inf]);
xlabel('time steps');