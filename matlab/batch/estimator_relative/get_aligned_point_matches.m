function cr=get_aligned_point_matches(options,idx)

ncams = length(idx);
cr = cell(ncams,1); % corners cell

disp('    ');
disp(['Processing images: ' num2str(idx(1)) ' to ' num2str(idx(ncams))]);

% get the images
im = get_bundle_images(options, idx);

% matching parameters (see elas.h for a full list)
% my_params;
% [D1,D2] = elasMex(im{1},im{2});%,param);
% D1
% imagesc(D1)
% pause
% fast corners in reference image only
%cprintf('blue','Fast corners\n');
fprintf('Fast corners\n');
[cr{1},use_r]=initialise_keyframe(im{1},options);
[cr,use_r,im]=track_bundle_points(cr,im,use_r,options); % track corners in the remaining images in the bundle
% also, minimum features disparity is applied here

%cprintf('blue','SIFT features\n');
fprintf('SIFT features\n');
[fr{1},d1,use_f]=initialise_keyframe_features(im{1},options);
[fr,use_f,im]=track_bundle_points_features(fr,d1,im,use_f,options);

% merge corners and features
use = [use_r; use_f];
for i=1:length(cr)
    cr{i}=[cr{i};fr{i}];
end

% plots
if 0
    for ii=2:length(cr)
        status=(cr{ii}(:,3)==1);
        clf; imshow([im{1},im{ii}]); hold on;
        plot(cr{1}(status,1),cr{1}(status,2),'+');
        plot(cr{ii}(status,1)+options.imgsize(2),cr{ii}(status,2),'+');
        pause
    end
end

% calibrate the points in all the images
% also, untracked or unwanted corners (use = 0) are removed here
cr = calibrate_bundle_points(cr, idx, use, options);

% plots
if 0
    colors = hsv(ncams);
    for i = 1:ncams
        status = cr{i}(3,:);
        plot(cr{i}(1,status==1), cr{i}(2,status==1), 'o', 'color', colors(i,:)); hold on;
        axis ij;
    end
    pause
end

%%%%%%%%
%%%%%%%%
%%%%%%%%
function im = get_bundle_images(options, idx)
ncams = length(idx);
im = cell(ncams,1); % images cell
%warning('left images were cropped');
%fprintf('images: ');
for i = idx
    %fprintf(['#',num2str(i), ', ']);
    if mod(i,2);
        I = imread(strcat(options.img_folder{1}, options.cam_left.image{(i+1)/2}));
        if size(I,3) == 1 % bayer decoding?
            I = demosaic(I, 'grbg');
        end
        %I = I(5:end-4,5:end-4,:);
    else
        I = imread(strcat(options.img_folder{2}, options.cam_right.image{(i+0)/2}));
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
%%%%%%%%
%%%%%%%%
%%%%%%%%
function [cr, use] = initialise_keyframe(im, options)
cr = fast9(im, options.fastthreshold, options.fastnonmax);
xlimits = [options.fastmargin,size(im,2)-options.fastmargin];
ylimits = [options.gridhorizon,size(im,1)-options.fastmargin];
cr = cr(cr(:,1)>xlimits(1)&cr(:,1)<xlimits(2)&...
        cr(:,2)>ylimits(1)&cr(:,2)<ylimits(2),:);
cr(:,3) = isfinite(cr(:,1))&isfinite(cr(:,2));
fprintf(['#' num2str(1)]);
fprintf(['(' num2str(sum(cr(:,3)==1)) '), ']);
%use=cr{1}(:,3);
use = zeros(size(cr,1),1); % added for more views triangulation fix
if 0
    clf; imshow(im); hold on;
    plot(cr(:,1),cr(:,2),'+');
    pause
end
%%%%%%%%
%%%%%%%%
%%%%%%%%
function [cr,use,im]=track_bundle_points(cr,im,use,options)
INITIALISE_NEXTPOINTS=0; % extract second image points and match before
% optical flow refinement
ncams=length(im);
% initialise previous points
prevPts=cell(1,size(cr{1},1));
for ii=1:size(cr{1},1) % put corners in optical flow cell point format
    prevPts{ii}=cr{1}(ii,1:2);
end
if INITIALISE_NEXTPOINTS==1
    %[~, desc1] = vl_sift(im{1},'frames',cr{1}(:,[1 2])); % Build a sift descriptor
    [desc1,~]=extractFeatures(im{1},cr{1}(:,[1 2]),'Method','SURF'); % Build a sift descriptor
    kdtree=vl_kdtreebuild(desc1','NumTrees',12); % Build kdtree (vl_sift)
end
%Criteria=struct('type','Count+EPS','maxCount',1000,'epsilon',eps);
Criteria=struct('type','EPS','maxCount',1e10,'epsilon',1e-10);
for ii=2:ncams
    fprintf(['#' num2str(ii)]);
    % initialise next points
    cr{ii}=cr{1}; % initialise to previous points
    if INITIALISE_NEXTPOINTS==1 % match features first usind descriptors?
        cr_ii=initialise_keyframe(im{ii},options);
        %[~, desc2]=vl_sift(im{ii},'frames',cr_ii); % Build a sift descriptor
        [desc2,~]=extractFeatures(im{ii},cr_ii(:,1:2),'Method','SURF'); % Build a descriptor
        matches=query_kdtree(kdtree,desc1',desc2');
        cr{ii}(matches(1,:),:)=cr_ii(matches(2,:),:); % update with matched next points
    end
    nextPts=cell(1,size(cr{ii},1));
    for jj=1:size(cr{ii},1) % put corners in optical flow cell point format
        nextPts{jj} = cr{ii}(jj,[1 2]);
    end
    % refine matching using optical flow?
    [nextPts,status,error]=cv.calcOpticalFlowPyrLK(im{1},im{ii},prevPts, ...
        'InitialFlow',      nextPts, ...
        'GetMinEigenvals',  true, ...
        'MinEigThreshold',  0.001, ... %.08
        'MaxLevel',         7, ...
        'WinSize',          options.LKWinSize, ...
        'Criteria',         Criteria);
    % cut the bundle if not enough corners were tracked
    %if sum(status==1&error<minerror) < options.mincorrnr;
    %    im = im(1:ii-1);
    %    cr = cr(1:ii-1);
    %    fprintf('\n');
    %    return;
    %end
    cr{ii}=vertcat(nextPts{:});
    % minimum features disparity and validity
    cr{ii}(:,3)=verify_point_track(status,cr{1},cr{ii},error,options);
    fprintf(['(' num2str(sum(cr{ii}(:,3)==1)) '), ']);
    if ~mod(ii, 10); fprintf('\n'); end;
    use = use | cr{ii}(:,3);
    % added for more views triangulation fix
    if 0
        status=(cr{ii}(:,3)==1);
        clf; imshow([im{1},im{ii}]); hold on;
        plot(cr{1}(status,1),cr{1}(status,2),'+');
        plot(cr{ii}(status,1)+options.imgsize(2),cr{ii}(status,2),'+');
        pause
    end
end
%%%%%%%%
%%%%%%%%
%%%%%%%%
function [cr,d1,use]=initialise_keyframe_features(im,options)
% compute keypoints
I=single(im);

[f1,d1] = vl_sift(I,'PeakThresh',0,'edgethresh',30, 'FirstOctave', -1);

%binSize = 8; magnif = .1;
%Is = vl_imsmooth(I, sqrt((binSize/magnif)^2 - .25)) ;
%[f1, d1] = vl_dsift(Is, 'size', binSize) ;

%[f1,d1]=vl_phow(im2single(im));

cr=zeros(size(f1,2),3);
cr(:,1:2)=f1(1:2,:)';
cr(:,3)=isfinite(cr(:,1))&isfinite(cr(:,2));
fprintf(['#' num2str(1)]);
fprintf(['(' num2str(sum(cr(:,3)==1)) '), ']);
%xlimits=[options.fastmargin size(im,2)-options.fastmargin];
%ylimits=[options.gridhorizon size(im,1)-options.fastmargin];
%idx=cr(:,1)>xlimits(1)&cr(:,1)<xlimits(2)&...
%    cr(:,2)>ylimits(1)&cr(:,2)<ylimits(2);
%cr=cr(idx,:);
%d1=d1(:,idx);
use=zeros(size(cr,1),1); % added for more views triangulation fix
if 0
    clf; imshow(im); hold on;
    plot(cr(:,1),cr(:,2),'+');
    pause
end
%%%%%%%%
%%%%%%%%
%%%%%%%%
function [cr,use,im]=track_bundle_points_features(cr,d1,im,use,options)
ncams=length(im);
cr{1}=cr{1}(:,1:3); % then assign only using location and status
for ii=2:ncams
    fprintf(['#' num2str(ii)]);
    cr{ii}=zeros(size(cr{1},1),3);
    
    I=single(im{ii});
    [f2,d2]=vl_sift(I,'PeakThresh',0,'edgethresh',30, 'FirstOctave', -1);
    %binSize = 8; magnif = 3;
    %Is = vl_imsmooth(I, sqrt((binSize/magnif)^2 - .25));
    %[f2,d2]=vl_dsift(Is,'size',binSize);
    
    %[f2,d2]=vl_phow(im2single(im{ii}));
    
    [matches,~]=vl_ubcmatch(d1,d2,1.5);
    cr{ii}(matches(1,:),1:2)=f2(1:2,matches(2,:))';
    % minimum features disparity and validity
    status=ones(size(cr{1},1),1);
    error=zeros(size(cr{1},1),1);
    cr{ii}(:,3)=verify_point_track(status,cr{1},cr{ii},error,options);
    fprintf(['(' num2str(sum(cr{ii}(:,3)==1)) '), ']);
    if ~mod(ii, 10); fprintf('\n'); end;
    use = use | cr{ii}(:,3);
    % added for more views triangulation fix
    if 0
        status=(cr{ii}(:,3)==1);
        clf; imshow([im{1},im{ii}]); hold on;
        plot(cr{1}(status,1),cr{1}(status,2),'+');
        plot(cr{ii}(status,1)+options.imgsize(2),cr{ii}(status,2),'+');
        pause
    end
end
%%%%%%%%
%%%%%%%%
%%%%%%%%
function status=verify_point_track(status,p1,p2,error,options)
% minimum features disparity and validity
minerror=1000; % very high, because is not considered at this stage
vis=remove_points_at_infinity(p1,p2,options.mindisp);
%use = use & (vis' & status & isfinite(cr{ii}(:,1)));
status=vis'...
    &status...
    &error<minerror...
    &isfinite(p2(:,1))&isfinite(p2(:,2))...
    &p2(:,1)<options.imgsize(2)&p2(:,2)<options.imgsize(1)...
    &p2(:,1)>0&p2(:,2)>0;
%%%%%%%%
%%%%%%%%
%%%%%%%%
function cr = calibrate_bundle_points(cr, idx, use, options)
ncams = length(idx);
for ii=1:ncams
    if size(cr{ii},2)>0
        cr{ii} = cr{ii}(use==1,:);  % remove points not visible in any of the images
        [K, kc] = get_intrinsics(options, idx(ii)); % remove lense distortion from the key-points
        x = remove_lens_distortion(cr{ii}(:,1:2), kc, K);
        x = K\pextend(x'); % calibrate
        cr{ii} = [x(1:2,:); cr{ii}(:,3)'];
    end
end
%%%%%%%%
%%%%%%%%
%%%%%%%%
function matches = query_kdtree(kdtree, desc1, desc2)
distRatio = 0.5; % less than 1, smaller means more strict (1/2 is more strict than 1/1.5)
index = vl_kdtreequery(kdtree, desc1, desc2, 'MAXCOMPARISONS', 50, 'NUMNEIGHBORS',2);
maxvals = sum(desc1(:,index(1,:)).*desc2, 1);
secondmaxvals = sum(desc1(:,index(2,:)).*desc2, 1);
matches = zeros(1,size(desc2,2));
matches(acos(maxvals) < distRatio*acos(secondmaxvals)) = ...
    index(1, acos(maxvals) < distRatio*acos(secondmaxvals));
ind2 = find(matches ~= 0);
ind1 = matches(ind2);
matches = [ind1; ind2];