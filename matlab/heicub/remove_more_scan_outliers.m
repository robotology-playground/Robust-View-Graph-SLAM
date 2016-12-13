function [scan, kpts] = remove_more_scan_outliers(options, pose, C, nviews, scan, kpts, ncam)
ncams = options.ncams;
nimages = size(scan, 2);
edges = vertcat(C.edge);
P = cell(1, size(pose, 2));
for i = 1:size(pose, 2)
    P{i} = [w2R(pose(4:6,i)), pose(1:3,i); 0 0 0 1];
end
for k = 1:nimages;

    if nargin<7  % this is because the old code didn't save ncam

        % check if this was done before for this camera
        if isempty(scan{k});
            continue;
        end
        % Adjust ncams for last few image frames
        if nimages-k+1 < ncams;
            ncams = nimages-k+1;
            options.ncams = ncams;
        end
        % Get conrners and tracks relative to the reference image
        p = get_aligned_point_matches(options, k);
        % Initialise constraints, and estimates
        [CC, CCt, sw, xs, sigmas, p] = initialise_inverse_depth(p, Pkin(k+(1:ncams)-1), options);
        % Get switching results
        i = find(edges(:,1)==k);
        sw = C(i(1)).sw;
        % Generate scans and remove outliers
        vis = zeros(1, size(p{1},2));
        cam = zeros(1, ncams);
        for i=1:length(CC)
            vis(CC(i).kpt) = vis(CC(i).kpt) + sw(i); %just add'm up ...
            cam(CC(i).cam) = cam(CC(i).cam) + sw(i);
        end
        vis = vis(vis>6);

    else

        vis = ncam{k}; % this is the new code ncam

    end

    index = find(vis > nviews);
    kpts{k} = kpts{k}(1:2, index);
    scan{k} = scan{k}(1, index);
    
    %clf; hold on;
    check_depth_map(scan(1:k), kpts(1:k), options, P(1:k), 1.5, 0, ones(1,k));
    axis([-2 .5 0 2]); grid on;
    drawnow;

end
