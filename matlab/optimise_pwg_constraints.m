function C = optimise_pwg_constraints(C,kpts,options)
%C=optimise_pwg_constraints(C,kpts,options)
%
% Optimises motion constraints in C using robust nonlinear least-squares
% This function uses the batch image optimisation class PwgOptimser
%
%	For all robots.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

% configurations
config_rswitch();

keyframe = 0;

%for k=1:length(C)
k = 0;
while true

	k = k + 1;
	
	% Check if end of constraints is reached, and hence optimise and exit.
	if k>length(C) || (is_new_keyframe(C(k),keyframe) && keyframe>0)
		ncams = length(camera);
		npts = length(tracks);
		xs = [xc;xf];
		assert(length(xs)==ncams*6+npts, 'state vector should be of length ncams*6+npts');
		if 0 % Debugging;
			plot_bundle_state(xs,Cimg,ncams);
			pause
		end
		x = optimise_image_data(xs,Cimg,sw,ncams,options);
		for i = 2:ncams
			C(idx(i-1)).t = x((i-1)*6+(1:3));
			C(idx(i-1)).a = x((i-1)*6+(4:6));
		end
		C(idx(1)).xf = x(ncams*6+(1:npts)); % FIXME: assumes cam1 only has a scan
		if k>length(C)
			break
		end
	end

	% Check if this is a new keyframe, and hence reset vectors
	if is_new_keyframe(C(k),keyframe) % a new edge to add
		keyframe = C(k).edge(1);
		xc = zeros(6,1);
		xf = [];
		camera = keyframe;
		tracks = C(k).matches(1,:);
		Cimg = []; % image constraints structure
		sw = []; % ON/OFF switching vector
		idx = []; % graph constraint index
		fprintf(['Keyframe ' num2str(keyframe) ', ']);
	end

	% Include a new camera and scan states
	[xc,camera] = push_camera_state(C(k),xc,camera);
	if has_scan(C(k))
		[xf,tracks] = push_scan_state(C(k),xf,tracks);
	end
	
	% Include new image measurements
	[p1,p2] = get_correspondence(C(k),kpts,options);
	n = length(Cimg);
	s = 0;
	for j = 1:size(p1,2) % Fields order is important (mex requirement)
		[member,position] = ismember(C(k).matches(1,j),tracks);
		if member==1
			s = s+1;
			Cimg(n+s).cam = length(camera); % camera number
			Cimg(n+s).kpt = position; % point track
			Cimg(n+s).p1 = p1(1:2,j); % projection in keyframe
			Cimg(n+s).z = p2(1:2,j); % projection in second camera
			Cimg(n+s).R = R./options.K1(1,1); % projection noise matrix (pixels/focal_length)
			sw(n+s) = 0; % FIXME: all constraints are initially set as OFF, but maybe not?
		end
	end
	
	% Count this constraint into the list
	idx = [idx k];
	
end %for k = 1:length(C) or while true
end %optimise_pwg_constraints()

% Additional functions
function flag = is_new_keyframe(C,keyframe)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	flag = logical( keyframe~=C.edge(1) );
end %is_new_keyframe()

function flag = has_scan(C)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	flag = logical( isfield(C,'xf') );
	if flag;
		flag = logical( ~isempty(C.xf) );
	end
end %has_scan()

function [x,c] = push_camera_state(C,x,c)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	x=[x;C.t;C.a];
	c=[c C.edge(2)];
end %push_camera_state()

function [x,t] = push_scan_state(C,x,t)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	% FIXME: scan merging, not implemented yet
	x=[x;C.xf];
	t=t;
end %push_scan_state()

function [p1,p2]=get_correspondence(C,kpts,options)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	p1=kpts{C.edge(1)}(C.matches(1,:),1:2);
	p1=calibrate_image_points(p1,options,C.edge(1))';
	p2=kpts{C.edge(2)}(C.matches(2,:),1:2);
	p2=calibrate_image_points(p2,options,C.edge(2))';
end %get_correspondence()

function p=calibrate_image_points(p,options,k)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
	[K,kc]=get_intrinsics(options,k); % remove lense distortion from the key-points
	x=remove_lens_distortion(p(:,1:2),kc,K);
	x(:,1)=(x(:,1)-K(1,3))/K(1,1);
	x(:,2)=(x(:,2)-K(2,3))/K(2,2);
	p(:,1:2)=x(:,1:2);
end %calibrate_image_points()

function x=optimise_image_data(xs,C,sw,ncams,options)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	% FIXME: future work, change PwgOptimiser to run as a Class with global parameters
	% configurations
	config_rswitch();
	if TRUST_INIT_LINKS
    	Ct = C(sw == 1);
    	C = C(sw == 0);
    	sw = sw(sw == 0);
	else
    	Ct = [];
	end
	npts=length(xs)-6*ncams;
	fprintf(['running optimisation using ' ... 
		num2str(ncams) ' cameras and ' num2str(npts) ' points.\n']);
	[x,~,~,sw] = PwgOptimiser(C,Ct,xs,sw,ncams,options);
	x=xs;
end %optimise_image_data()

function plot_bundle_state(x,C,ncams)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
	plot_pose_graph_w(x(1:ncams*6));
	idx=0;
	for k=1:length(C)
		if C(k).cam==2 % FIXME: assumes the second camera is used to construct the map
			idx=idx+1;
			p(:,idx)=C(k).p1;
			r(1,idx)=1/x(ncams*6 + C(k).kpt);
		end
	end
	get_scan_from_range(p,r,1);
end %plot_bundle_state()

function plot_pose_graph_w(x, i)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2015
	cols = 'krbgcmy';
	persistent IDX
	if isempty(IDX) || IDX == length(cols), IDX = 0; end
	IDX = IDX + 1;
	if nargin == 2, IDX = mod(i-1, length(cols)) + 1; end

	plot3(x(1:6:end), x(3:6:end), x(2:6:end), cols(IDX))
	hold on

	frame = .001*[-0.1 1 nan,    0 0 nan, nan    0 0;
            0 0 nan, -0.2 1 nan, nan    0 0;
            0 0 nan,    0 0 nan, nan -0.2 1];
	idx = 1:6;
	while idx(1) < length(x)
    	fg = transform_to_global_w(frame, x(idx));
    	idx = idx + 6;
    	plot3(fg(1,:),fg(3,:),fg(2,:), cols(IDX))
	end
	axis equal, grid on
	drawnow
end %plot_pose_graph_w()

function p = get_scan_from_range(p, r, verbose)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2015
	if size(p,1) < 3; p = pextend(p); end;
	rim = sqrt(p(1,:).*p(1,:) + p(2,:).*p(2,:) + p(3,:).*p(3,:));
	d = r./rim;
	p(1,:) = d.*p(1,:);
	p(2,:) = d.*p(2,:);
	p(3,:) = d.*p(3,:);
	if verbose % Debugging; 
		plot3(p(1,:),p(2,:),p(3,:),'+');
		axis equal, grid on
		drawnow
	end
end %get_scan_from_range()
