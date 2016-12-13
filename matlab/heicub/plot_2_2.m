function plot_2(scan, ncams, npts, C, Pkin, kpts, options)
use = zeros(1, size(scan,2));

% idx = [7 9 21];
% use(idx) = 1;

for tt = 1 : 2 : size(use,2)
    use(tt) = sum(npts{tt}>10) > (.75*ncams); % this avoids using
    % scans from bad bundles, for example, when the robot isnt
    % moving. You probably need to not process these bundles to start
    % with, rather than spending time on them. However, this is a
    % temporary hack.
end

disp([num2str(sum(use)), ' scans has been generated']);
if sum(use)
    clf; hold on;
    
    % plot the scans
    check_depth_map(scan, kpts, options, Pkin, 5, 0, use, 2);
    
    % plot trajectories estimated (from C)
    edges = vertcat(C.edge);
    refs = unique(edges(:,1)); % key-frames
    for ii = 1 : length(refs)
        idx = find(edges(:,1)==refs(ii)); % images wrt a reference frame
        pose = zeros(6,length(idx)+1); % +1 for the reference frame
        cam = ones(1,length(idx)+1);
        cam(1) = mod(refs(ii),2);
        for jj = 1 : length(idx)
            pose(:,jj+1) = C(idx(jj)).z; % +1 for the reference frame
            cam(jj+1) = mod(C(idx(jj)).edge(2),2); % flag, left or right ?
        end
        ref_pose = camera_matrix_to_pose(Pkin{refs(ii)}); % reference frame
        pose = transform_to_global_w(pose,ref_pose); % other frames
        %plot3(pose(2,:),pose(1,:),pose(3,:),'r*','makersize',2);
        plot3(pose(2,cam==1),pose(1,cam==1),pose(3,cam==1),...
            'r','linewidth',1); % left camera trajectory
        plot3(pose(2,cam==0),pose(1,cam==0),pose(3,cam==0),...
            'g','linewidth',1); % right camera trajectory
    end
    
    % plot trajectories from kinematics
    %pose = camera_matrix_to_pose(Pkin(1:max(edges(:)))); % segment only ?
    pose = camera_matrix_to_pose(Pkin); % all kinematics trajectory ?
    plot3(pose(2,1:2:end),pose(1,1:2:end),pose(3,1:2:end), ...
        'color',[.75 .75 .75],'linewidth',1); % left kinematics
    plot3(pose(2,2:2:end),pose(1,2:2:end),pose(3,2:2:end), ...
        'color',[.75 .75 .75],'linewidth',1); % right kinematics
    
    axis equal; axis tight; grid on; box on;
    xlabel('Eastern (m)'); ylabel('Forward (m)');
    drawnow;
end

function x = camera_matrix_to_pose(P)
if iscell(P);
    n = length(P);
    x = zeros(6,n);
    for i = 1 : n
        x(:,i) = [P{i}(1:3,4);R2w(P{i}(1:3,1:3))'];
    end
else
    x = [P(1:3,4); R2w(P(1:3,1:3))'];
end
