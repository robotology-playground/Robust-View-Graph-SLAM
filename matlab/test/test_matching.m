

% match features
% Tariq Abuhashim - August 2014, iCub

function [matches, mpts1, mpts2] = test_matching(kpts1, kpts2, desc1, desc2, options)

%fprintf('matching .... ');
t_start = tic;

% Build kdtree (kaze)
if strcmp(options.method, 'kaze');
    
    % Build kdtree (vl_sift)
    kdtree = vl_kdtreebuild(desc1,'NumTrees', 12);
    distRatio = options.kazeratio; % less than 1, smaller means more strict (1/2 is more strict than 1/1.5)
    index = vl_kdtreequery(kdtree, desc1, desc2, 'MAXCOMPARISONS', 50, 'NUMNEIGHBORS',2);
    maxvals = sum(desc1(:,index(1,:)).*desc2, 1);
    secondmaxvals = sum(desc1(:,index(2,:)).*desc2, 1);
    matches = zeros(1,size(desc2,2));
    matches(acos(maxvals) < distRatio*acos(secondmaxvals)) = ...
        index(1, acos(maxvals) < distRatio*acos(secondmaxvals));
    ind2 = find(matches ~= 0);
    ind1 = matches(ind2);
    %[~,I,~]=unique(ind1); % remove non uniqe matches (Caused by NN approximation?)
    %ind2 = ind2(I); ind1 = ind1(I);
    matches = [ind1; ind2];
    mpts1 = kpts1(matches(1,:),:)'; % 2 x N
    mpts2 = kpts2(matches(2,:),:)'; % 2 x N
    
    % % vl_ubcmatch
    % %THRESH = 10; % 3 with stereo setting, 10 for monocular
    % [matches, ~] = vl_ubcmatch(desc1, desc2, THRESH);
    % % A descriptor desc1 is matched to a descriptor desc2 only if the distance
    % % d(desc1,desc2) multiplied by THRESH is not greater than the distance of
    % % desc1 to all other descriptors. The default value of THRESH is 1.5
    % mpts1 = kpts1(matches(1,:), :)'; % 2 x N
    % mpts2 = kpts2(matches(2,:), :)'; % 2 x N
    
elseif strcmp(options.method, 'sift');
    
    % % Build kdtree (vl_sift)
    % kdtree = vl_kdtreebuild(desc1,'NumTrees', 12);
    % [index, dist] = vl_kdtreequery(kdtree, desc1, desc2, 'MAXCOMPARISONS', 50,...
    %     'NUMNEIGHBORS',1);
    % ind2 = find(index ~= 0);
    % ind1 = index(ind2);
    % [~,I,~]=unique(ind1);
    % ind2 = ind2(I);
    % ind1 = ind1(I);
    % matches = [ind1; ind2];
    % mpts1 = kpts1(matches(1,:), :)'; % 2 x N
    % mpts2 = kpts2(matches(2,:), :)'; % 2 x N
    
    % vl_ubcmatch
    THRESH = options.siftratio; % 3 with stereo setting, 10 for monocular
    [matches, ~] = vl_ubcmatch(desc1, desc2, THRESH);
    % A descriptor desc1 is matched to a descriptor desc2 only if the distance
    % d(desc1,desc2) multiplied by THRESH is not greater than the distance of
    % desc1 to all other descriptors. The default value of THRESH is 1.5
    mpts1 = kpts1(matches(1,:), :)'; % 2 x N
    mpts2 = kpts2(matches(2,:), :)'; % 2 x N
    
end

fprintf(['matches ',num2str(size(matches, 2))]);
t_end = toc(t_start); 
fprintf([',  ', num2str(t_end), 's\n']);