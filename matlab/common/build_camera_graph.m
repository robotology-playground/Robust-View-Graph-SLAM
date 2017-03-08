function [C,kpts]=build_camera_graph(options)
%[C,kpts]=build_camera_graph(options)
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

nimages=size(options.cam_left.image,2)+size(options.cam_right.image,2);
kpts=cell(1,nimages); desc=cell(1,nimages); C=[]; k=0;
thickness=options.ncams-1;
%nkeys=0;
ncorners=options.mincorrnr;
mindisp=options.mindisp;

fprintf('Building the camera graph ');
if isfield(options,'kazethreshold')
	fprintf(['(detector = KAZE<',num2str(options.kazethreshold),',',num2str(options.kazeratio),'>']);
end
fprintf([', thickness = ',num2str(thickness)]);
fprintf([', mincorrnr = ',num2str(ncorners)]);
fprintf([', mindisp = ',num2str(mindisp)]);
fprintf([') ...\n']);
for i=1:2:nimages-1	% FIXME: No edge from even to odd
					% This was assumed because even-to-odd edge results 
					% were strange
	%nkeys=nkeys+1;
	%if nkeys>options.nkeys
	%	break
	%end
    if isempty(kpts{i})
        [kpts{i},desc{i}]=get_akaze(options,i);
    end
    for j=i+1:min(i+thickness,nimages)
        if isempty(kpts{j})
            [kpts{j},desc{j}]=get_akaze(options,j);
        end
        matches=vl_ubcmatch(desc{i}',desc{j}',options.kazeratio); % 2.5
        vis=remove_points_at_infinity(kpts{i}(matches(1,:),1:2)', kpts{j}(matches(2,:),1:2)', mindisp);
        matches=matches(1:2,vis==1);
        %if size(matches,2)>ncorners
            k=k+1;
            C(k).edge=[i,j];	% FIXME: doesn't account for uncomputed edges
            C(k).weight=size(matches,2);
            C(k).matches=matches;
	    	fprintf('Edge [%d,%d,%d].\n',i,j,size(matches,2));
        %end
    end
end
if isfield(options,'save');
   save(strcat(options.save,'/camera_graph.mat'),'-v7.3','C','options','kpts');
end
end
%
%
% Additional functions
function [c,d]=get_akaze(options,j)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
	im=get_image_k(options,j);
	detector=cv.FeatureDetector('KAZE','Threshold',options.kazethreshold);
	keypoints=detector.detect(im);
	c=vertcat(keypoints.pt);
	c(:,3)=c(:,1)<options.imgsize(2)&c(:,2)<options.imgsize(1)&c(:,1)>0&c(:,2)>0;
	extractor=cv.DescriptorExtractor('SIFT');
	d=extractor.compute(im,keypoints);
end

function im=get_image_k(options,k)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
	if mod(k,2)==1
    	I=imread(strcat(options.img_folder{1},options.cam_left.image{(k+1)/2}));
    	if options.splitimage==1 % split the stereo image (for R1 composite images)
        	[m,n,k]=size(I);
        	I=I(1:m,1:n/2,1:k);
    	end
    	if size(I,3)==1 % bayer decoding?
        	I=demosaic(I,'grbg');
    	end
	else
    	I=imread(strcat(options.img_folder{2},options.cam_right.image{(k+0)/2}));
    	if options.splitimage==1 % split the stereo image (for R1 composite images)
    	    [m,n,k]=size(I);
        	I=I(1:m,n/2+(1:n/2),1:k);
    	end
    	if size(I,3)==1 % bayer decoding?
        	I=demosaic(I,'grbg');
    	end
	end
	if length(size(I))==3
	    %I=I(:,:,2);
	    I=rgb2gray(I);
	end
	im=I;
end
