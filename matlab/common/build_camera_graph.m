function C=build_camera_graph(options)
%C=build_camera_graph(options)
%
% Builds a graph with a given thickness and populates edges with 
% weights using frame correspondences.
%
%	For all robots.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

nimages = size(options.cam_left.image,2)+size(options.cam_right.image,2);
kpts = cell(1,nimages); 
desc = cell(1,nimages);
thickness = options.ncams-1;
ncorners = options.mincorrnr;
mindisp = options.mindisp;
C = []; k = 0; %nkeys=0;
 
fprintf(' - Building the camera graph ');
fprintf(['(detector = ',options.detector,'<',strjoin(arrayfun(@(x) num2str(x),options.detector_param,'UniformOutput',false),','),'>']);
fprintf([', thickness = ',num2str(thickness)]);
fprintf([', mincorrnr = ',num2str(ncorners)]);
fprintf([', mindisp = ',num2str(mindisp)]);
fprintf([') ...\n']);

for i=1:2:nimages-1	% FIXME: No edge from even to odd
					% This was assumed because even-to-odd
					% edge results were strange
					
	%nkeys=nkeys+1;
	%if nkeys>options.nkeys
	%	break
	%end
    
    % process the reference frame
    switch options.detector
    case 'KAZE'
    	if isempty(kpts{i})
        	[kpts{i},desc{i}] = get_akaze(options,i);
        end
        c1 = kpts{i}; d1 = desc{i}; % KAZE is unique, so no need to repeat
	case 'FAST'
		c1 = get_fast(options,i); % Corners in each image are unique, but not their tracks
	end % switch

	% track in the remaining frames
    for j=i+1:min(i+thickness,nimages)
    
    	switch options.detector
    	case 'KAZE'
        	if isempty(kpts{j})
            	[kpts{j},desc{j}] = get_akaze(options,j);
        	end
        	c2 = kpts{j}; d2 = desc{j}; % KAZE is unique, so no need to repeat
        	kazeratio = options.detector_param(2);
        	matches = vl_ubcmatch(desc{i}',desc{j}',kazeratio); % Matching using vlFeat
        case 'FAST'
        	c2 = track_fast(options,c1,i,j);
        	% FIXME: fast and optical flow are not unique, we need to find away to merge those in images
        	% i.e., a corner extracted from frame 2 should be merged with the same point tracked into frame 2
        	% for example, using nearest-neighbour, or buy extracting local descriptors.
        	% Hence, c2 and kpts{j} should have the same ID
        	matches = [1:size(c1,1);1:size(c2,1)]; % tracked in the same order using optical flow
        end % switch
        
        vis = remove_points_at_infinity(c1(matches(1,:),1:2)',c2(matches(2,:),1:2)',mindisp);
        matches = matches(1:2,vis==1);
        
        if size(matches,2)>ncorners
            k = k+1;
            C(k).edge = [i,j]; % FIXME: doesn't account for uncomputed edges
            C(k).weight = size(matches,2);
            C(k).kpts1 = c1(matches(1,:),:);
            C(k).kpts2 = c2(matches(2,:),:);
            C(k).matches = matches;
	    	fprintf('Edge [%d,%d,%d].\n',i,j,size(matches,2));
	    	if 0
	    		clf;
	    		im1 = get_image_k(options,i);
				im2 = get_image_k(options,j);
				imshow([im1,im2]); hold on;
				plot(c1(matches(1,:),1),c1(matches(1,:),2),'r+');
				plot(c2(matches(2,:),1)+size(im1,2),c2(matches(2,:),2),'b+');
	   			drawnow;
	    	end
        end % if
        
    end % for
    
end
if isfield(options,'save');
   save(strcat(options.save,'/camera_graph.mat'),'-v7.3','C','options','kpts');
end
end
%
%
% Additional functions
%
%
function [c,d]=get_akaze(options,j)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
	threshold=options.detector_param(1);
	im=get_image_k(options,j);
	detector=cv.FeatureDetector('KAZE','Threshold',threshold);
	keypoints=detector.detect(im);
	c=vertcat(keypoints.pt);
	c(:,3)=c(:,1)<options.imgsize(2)&c(:,1)>0&... % status of each point feature
			c(:,2)<options.imgsize(1)&c(:,2)>0; % 1-valid, 0-invalid
	extractor=cv.DescriptorExtractor('SIFT');
	d=extractor.compute(im,keypoints);
end
%
%
function cr=get_fast(options,k)
	margin=options.detector_param(1);
	threshold=options.detector_param(2);
	nonmax=options.detector_param(3);
	horizon=options.gridhorizon;
	im=get_image_k(options,k);
	%cr=fast9(im,threshold,nonmax);
	detector=cv.FeatureDetector('FastFeatureDetector', ...
		'Threshold',threshold,'NonmaxSuppression',nonmax);
	keypoints=detector.detect(im);
	cr=vertcat(keypoints.pt);
	xlimits=[margin,size(im,2)-margin];
	ylimits=[horizon,size(im,1)-margin];
	cr=cr(cr(:,1)>xlimits(1)&cr(:,1)<xlimits(2)&...
      	cr(:,2)>ylimits(1)&cr(:,2)<ylimits(2),:);
	cr(:,3)=isfinite(cr(:,1))&isfinite(cr(:,2));
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
	winsize=options.detector_param(4);
	Criteria=struct('type','EPS','maxCount',1e10,'epsilon',1e-10) ;
	[nextPts,status,error]=cv.calcOpticalFlowPyrLK(im1,im2,prevPts, ...
		'GetMinEigenvals',  true, ...
		'MinEigThreshold',  0.001, ... %.08
		'MaxLevel',         7, ...
		'WinSize',          [winsize, winsize], ...
		'Criteria',         Criteria) ;
	% filter tracks (status, minimum features disparity and validity)
	c2=vertcat(nextPts{:});
	c2(:,3)=verify_point_track(c1,c2,options,status,error);
end
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
end
%
%
function im=get_image_k(options,k)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
%
% this functions reads an image given left and right image
% names. The function matches the way names are
% saved in options, which can be adapted differently if needed.

	if mod(k,2)==1
    	I=imread(strcat(options.img_folder{1},options.cam_left.image{(k+1)/2}));
    	if options.splitimage==1 % split the stereo image (for composite images)
        	[m,n,k]=size(I);
        	I=I(1:m,1:n/2,1:k);
    	end
    	if size(I,3)==1 % bayer decoding?
        	I=demosaic(I,'grbg'); % iCub implements grbg coding, change to your case
    	end
	else
    	I=imread(strcat(options.img_folder{2},options.cam_right.image{(k+0)/2}));
    	if options.splitimage==1 % split the stereo image (for composite images)
    	    [m,n,k]=size(I);
        	I=I(1:m,n/2+(1:n/2),1:k);
    	end
    	if size(I,3)==1 % bayer decoding?
        	I=demosaic(I,'grbg'); % iCub implements grbg coding, change to your case
    	end
	end
	
	if length(size(I))==3
	    %I=I(:,:,2); % use green channel only ?
	    I=rgb2gray(I); % convert to gray-scale ?
	end
	im=I;
	
	% one of the datasets has left images with 488x648 resolution.
	% remove this if statement otherwise
	if size(im,1)==488&size(im,2)==648
		im=im(5:end-4,5:end-4);
	end
	
end
