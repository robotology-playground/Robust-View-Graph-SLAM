function check_depth_map(scan,kpts,options,P,MAX_RANGE,calib,use,root)

idx = [];
for i = 1:length(scan)
    if ~isempty(scan{i}) && use(i) == 1
        idx = [idx i];
    end
end
color(idx,:) = hsv(length(idx));
%i = randperm(count);
%color = color(i,:);

if nargin < 8; root = 0; end;

for i = idx
    %if ~isempty(scan{i}) && use(i)
        
        xf = scan{i};
        p = kpts{i}(1:2,:);
        
        % calibrate the points ?
        if calib == 1;  
            [K, k] = get_intrinsics(options, i);
            p = remove_lens_distortion(p', k, K);
            p = K\pextend(p');
        end
        
        % transform from depth/inverse depth to 3d scans
        idx = xf > 0; % remove negative inverse depth
        r = 1./xf(idx);
        X = get_scan_from_range(p(1:2,idx), r);
        X = X(:, r<MAX_RANGE); % remove distant points
        
        % transform to global camera coordinates
        pose = camera_matrix_to_pose(P{i});
        X = transform_to_global_w(X, pose);
        
        % plotting
        hold on;
        if root == 1
            plot3(X(1,:), X(3,:), -X(2,:), '.', 'color', ...
                color(i,:),'markersize',2); % left camera coordinates
        elseif root == 2
            plot3(-X(2,:), X(1,:), X(3,:), '.', 'color', ...
                color(i,:),'markersize',2); % floating-base coordinates
        else
            plot3(-X(3,:), X(2,:), X(1,:), '.', 'color', ...
                color(i,:),'markersize',2); % waist coordinates
        end
        axis equal; grid on; view(2); %drawnow;

    %end
end
%
%
function p = get_scan_from_range(p, r)
rim = sqrt(sum(p(1:2,:).^2) + 1);
d = r./rim;
p(1,:) = p(1,:).*d;
p(2,:) = p(2,:).*d;
p(3,:) = d;
%
%
function x = camera_matrix_to_pose(P)
if iscell(P);
    n=length(P);
    x=zeros(6,n);
    for i=1:n
        x(:,i)=[P{i}(1:3,4);R2w(P{i}(1:3,1:3))'];
    end
else
    x=[P(1:3,4);R2w(P(1:3,1:3))'];
end