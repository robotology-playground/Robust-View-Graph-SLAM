function cr = get_aligned_point_matches_v2(options, ref_frame, no_cams)

% k : reference frame
% j : second frame

cr = struct('c',[],'x',[],'y',[],'s',[]); % corners cell
% c : camera id
% x,y : pixel coordinates
% s : status
nimages = 2*length(options.cam_left.image); % total available number of images
j = ref_frame; % global index of this reference frame
k = 1; % local index of this reference frame

disp('    ');
cprintf('Green',['Processing key-frame: ' num2str(j) '\n']);

% find the corners
%c1=get_fast(options,j);
%[c1,d1]=get_sift(options,j);
[c1,d1] = get_akaze(options,j); % status set using image limits
c1_ = calibrate_image_points(c1,options,j); % no change on status
cr(k).c = j;
cr(k).x = c1_(:,1)';
cr(k).y = c1_(:,2)';
cr(k).s = c1_(:,3)'; % same status as in c1
vis = cr(k).s; % count number of cameras seeing each point
if options.verbose>0
    fprintf(['#' num2str(j) '(' num2str(sum(cr(k).s==1)) '), ']);
end

% find the tracks
%while k<ncams && j<nimages % get at least ncams frames
while (j-ref_frame+1)<no_cams && j<nimages % only look at the next no_cams frames
    j = j+1;
    %c2 = track_fast(options,c1,ref,j);
    %c2 = track_sift(options,d1,j);
    c2 = track_akaze(c1,d1,options,j); % status is set based on visibility
    c2_ = calibrate_image_points(c2,options,j); % no change on status
    if sum(c2_(:,3))>options.mincorrnr
        k = k+1; % add a new track
        cr(k).c = j; % local index of every frame
        cr(k).x = c2_(:,1)';
        cr(k).y = c2_(:,2)';
        cr(k).s = c2_(:,3)'; % same status as in c2
        vis = vis+cr(k).s; % count number of cameras seeing each point
        if options.verbose>0
            fprintf(['#' num2str(j) '(' num2str(sum(cr(k).s==1)) '), ']);
            if ~mod(k,10);
                fprintf('\n');
            end; % fill display by 10 results
        end
        if 0
            status = cr(1).s==1 & cr(k).s==1;
            clf; plot(cr(1).x(status),cr(1).y(status),'o');
            hold on; plot(cr(k).x(status),cr(k).y(status),'ro');
            title(num2str(j)); pause;
        end
    end
end
fprintf('\n');

% resolve tracks
nviews = 2; % minimum number of image point projections in the sequence.
for k = 1 : length(cr)
    cr(k).x=cr(k).x(vis>=nviews);
    cr(k).y=cr(k).y(vis>=nviews);
    cr(k).s=cr(k).s(vis>=nviews);
end

% Check data associations between frames 1 and 2
if 0
    
    i=1; j=3;
    
    im1 = get_image_k(options,i);
    im2 = get_image_k(options,j);
    clf; imshow([im1 im2]); hold on;
    
    status = cr(i).s==1 & cr(j).s==1;
    
    x = cr(i).x(status);
    y = cr(i).y(status);
    [K,kc] = get_intrinsics(options,i);
    [x,y] = UndistPointFor(x, y, kc(1), kc(2), kc(3), kc(4), kc(5));
    x = K(1,1)*x + K(1,3);
    y = K(2,2)*y + K(2,3);
    p1 = [x;y];
    plot(p1(1,:),p1(2,:),'go');
    
    x = cr(j).x(status);
    y = cr(j).y(status);
    [K,kc] = get_intrinsics(options,j);
    [x,y] = UndistPointFor(x, y, kc(1), kc(2), kc(3), kc(4), kc(5));
    x = K(1,1)*x + K(1,3);
    y = K(2,2)*y + K(2,3);
    p2 = [x+size(im1,2);y];
    plot(p2(1,:),p2(2,:),'go');
    
    line([p1(1,:);p2(1,:)],[p1(2,:);p2(2,:)],'color',[0 1 0]);
    axis tight
    
    pause
    
end

if 0
    ncams = length(cr);
    colors = hsv(ncams);
    for i = 1:ncams
        plot(cr(k).x(cr(k).s==1), cr(k).y(cr(k).s==1), 'o', 'color', colors(i,:)); hold on;
        axis ij;
    end
    pause
end
%
%
function cr=get_fast(options,k)
im=get_image_k(options,k);
cr=fast9(im,options.fastthreshold,options.fastnonmax);
xlimits=[options.fastmargin,size(im,2)-options.fastmargin];
ylimits=[options.gridhorizon,size(im,1)-options.fastmargin];
cr=cr(cr(:,1)>xlimits(1)&cr(:,1)<xlimits(2)&...
      cr(:,2)>ylimits(1)&cr(:,2)<ylimits(2),:);
cr(:,3)=isfinite(cr(:,1))&isfinite(cr(:,2));
if 0
    clf; imshow(im); hold on;
    plot(cr(:,1),cr(:,2),'+');
    pause
end
%
%
function c2=track_fast(options,c1,i,j)
im1=get_image_k(options,i);
im2=get_image_k(options,j);
% initialise previous points (in openCV format)
prevPts=cell(1,size(c1,1));
for ii = 1 : size(c1,1) % put corners in optical flow cell point format
    prevPts{ii}=c1(ii,1:2);
end
% tracking using optical flow?
Criteria=struct('type','EPS','maxCount',1e10,'epsilon',1e-10) ;
[nextPts,status,error]=cv.calcOpticalFlowPyrLK(im1,im2,prevPts, ...
    'GetMinEigenvals',  true, ...
    'MinEigThreshold',  0.001, ... %.08
    'MaxLevel',         7, ...
    'WinSize',          options.LKWinSize, ...
    'Criteria',         Criteria) ;
% filter tracks (status, minimum features disparity and validity)
c2=vertcat(nextPts{:});
c2(:,3)=verify_point_track(c1,c2,options,status,error);
if 0
    status=(cr{ii}(:,3)==1);
    clf; imshow([im{1},im{ii}]); hold on;
    plot(cr{1}(status,1),cr{1}(status,2),'+');
    plot(cr{ii}(status,1)+options.imgsize(2),cr{ii}(status,2),'+');
    pause
end
%
%
function [c,d]=get_sift(options,j)
im=get_image_k(options,j);
[f,d]=vl_sift(im2single(im));
c=f(1:2,:)';
c(:,3)=1;
d=double(d);
%
%
function c2=track_sift(options,d1,j)
[f2,d2]=get_sift(options,j);
matches=vl_ubcmatch(d1,d2,3);
c2=zeros(size(d1,1),3);
c2(matches(1,:),1:3)=f2(matches(2,:),1:3);
%
%
function [c,d]=get_akaze(options,j)
% git clone https://github.com/pablofdezalc/akaze.git
% hacks to compile akaze
% 1- in the file ../mex/akaze.ccp, add on top 
%       using namespace libAKAZE;
% 2- cv::imwrite causing undefined reference error at lines 1370 and 1388
%    in the file ../src/lib/AKAZE.cpp, comment out both lines
im=get_image_k(options,j);
%[c,d]=akaze(im,'dthreshold',.0001,'descriptor',3);
detector=cv.FeatureDetector('KAZE','Threshold',.0001); %.0001
%detector=cv.FeatureDetector('BRISK','Threshold',30,'PatternScale',1);
keypoints=detector.detect(im);
c=vertcat(keypoints.pt);
c(:,3)=c(:,1)<options.imgsize(2)&c(:,2)<options.imgsize(1)&c(:,1)>0&c(:,2)>0;
extractor=cv.DescriptorExtractor('SIFT');
%extractor=cv.DescriptorExtractor('SURF');
%extractor=cv.DescriptorExtractor('BRISK');
d=extractor.compute(im,keypoints);
%
%
function c2=track_akaze(c1,d1,options,j)
[f2,d2]=get_akaze(options,j);
%kdtree=vl_kdtreebuild(d1,'NumTrees', 12);
%distRatio=.65; % less than 1, smaller means more strict (1/2 is more strict than 1/1.5)
%index=vl_kdtreequery(kdtree,d1,d2,'MAXCOMPARISONS',50,'NUMNEIGHBORS',2);
%maxvals=sum(d1(:,index(1,:)).*d2,1);
%secondmaxvals=sum(d1(:,index(2,:)).*d2,1);
%matches=zeros(1,size(d2,2));
%matches(acos(maxvals)<distRatio*acos(secondmaxvals))=...
%index(1,acos(maxvals)<distRatio*acos(secondmaxvals));
%ind2=find(matches ~= 0);
%ind1=matches(ind2);

matches=vl_ubcmatch(d1',d2',2.5); %2.5
%c2=zeros(size(d1,1),3);
c2(:,1:2)=c1(:,1:2); c2(:,3)=0; 
c2(matches(1,:),1:3)=f2(matches(2,:),1:3);

% matcher=cv.DescriptorMatcher('FlannBasedMatcher', ...
%     'Index',{'KDTree','Trees',12},'Search',{'Sorted', true});
% %matcher=cv.DescriptorMatcher('BFMatcher','NormType','L2');
% matcher.add(d2);
% matcher.train(); % Optional for BruteForce matcher
% matches=matcher.match(d1);
% ind1=vertcat(matches.queryIdx)+1; % starts from 0
% ind2=vertcat(matches.trainIdx)+1; % starts from 0
% dist=vertcat(matches.distance);
% c2=zeros(size(d1,1),3);
% %c2(ind1(dist<300),1:3)=f2(ind2(dist<300),1:3);
% c2(ind1,1:3)=f2(ind2,1:3);
c2(:,3)=verify_point_track(c1(:,1:2),c2(:,1:2),options,c2(:,3));
%
%
function im=get_image_k(options,k)
if mod(k,2);
    I=imread(strcat(options.img_folder{1},options.cam_left.image{(k+1)/2}));
    if options.splitimage % split the stereo image
        [m,n,k]=size(I);
        I=I(1:m,1:n/2,1:k);
    end
    if size(I,3)==1 % bayer decoding?
        I=demosaic(I,'grbg');
    end
    %I=I(5:end-4,5:end-4,:);
else
    I=imread(strcat(options.img_folder{2},options.cam_right.image{(k+0)/2}));
    if options.splitimage % split the stereo image
        [m,n,k]=size(I);
        I=I(1:m,n/2+(1:n/2),1:k);
    end
    if size(I,3)==1 % bayer decoding?
        I=demosaic(I,'grbg');
    end
end
if length(size(I))==3
    %I = I(:,:,2);
    I=rgb2gray(I);
end
im=I;
%
%
function status=verify_point_track(p1,p2,options,tracked,error)
% minimum features disparity and validity
minerror=1000; % very high, because is not considered at this stage
vis=remove_points_at_infinity(p1,p2,options.mindisp);
status=vis' ...
    &isfinite(p2(:,1))&isfinite(p2(:,2)) ...
    &p2(:,1)<options.imgsize(2)&p2(:,2)<options.imgsize(1) ...
    &p2(:,1)>0&p2(:,2)>0;
if nargin>3
    status=status&tracked;
end
if nargin>4
    status=status&error<minerror;
end
%
%
function cr=calibrate_image_points(cr,options,k)
use=find(cr(:,3)==1); % only calibrate poists with status=1 (saves time)
if ~isempty(use)
    [K,kc]=get_intrinsics(options,k); % remove lense distortion from the key-points
    x=remove_lens_distortion(cr(use,1:2),kc,K);
    x(:,1)=(x(:,1)-K(1,3))/K(1,1);
    x(:,2)=(x(:,2)-K(2,3))/K(2,2);
    %x=K\pextend(x'); % calibrate
    cr(use,1:2)=x(:,1:2);
end