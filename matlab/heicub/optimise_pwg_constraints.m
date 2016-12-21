function C = optimise_pwg_constraints(C,kpts,options)
%C=optimise_pwg_constraints(C,kpts,options)
%
% Optimises motion constraints in C using robust nonlinear least-squares
% This function uses the batch image optimisation class PwgOptimser
%
% Tariq Abuhashim - 2016
% t.abuhashim@gmail.com
%
% iCub - Koroibot

% configurations
config_rswitch();

keyframe = 0;

for k = 1:length(C)

	if is_new_keyframe(C(k),keyframe) % a new edge to add

		if keyframe % FIXME: Last bundle is not optimised
			ncams = length(camera);
			npts = length(tracks);
			fprintf(['running optimisation using ' ... 
					num2str(ncams) ' cameras and ' num2str(npts) ' points.\n']);
		end

		keyframe = C(k).edge(1);
		xc = zeros(6,1);
		camera = keyframe;
		xf = [];
		tracks = C(k).matches(1,:);
		Cimg = [];
		fprintf(['Keyframe ' num2str(keyframe) ', ']);

	end

	[xc,camera] = push_camera_state(C(k),xc,camera);
	if has_scan(C(k))
		[xf,tracks] = push_scan_state(C(k),xf,tracks);
	end
	
	[p1,p2] = get_correspondence(C(k),kpts,options);
	n = length(Cimg);
	for j = 1:size(p1,2) % Fields order is important (mex requirement)
		Cimg(n+j).cam = length(camera); % camera number
		[member,position] = ismember(C(k).matches(1,j),tracks);
		Cimg(n+j).kpt = position; % point track
		Cimg(n+j).p1 = p1(1:2,j); % projection in keyframe
		Cimg(n+j).z = p2(1:2,j); % projection in second camera
		Cimg(n+j).R = R./options.K1(1,1); % projection noise matrix (pixels/focal_length)
	end

end %for k = 1:length(C)
end %optimise_pwg_constraints()

function flag = is_new_keyframe(C,keyframe)
	flag = logical( keyframe~=C.edge(1) );
end %is_new_keyframe()

function flag = has_scan(C)
	flag = logical( isfield(C,'xf') );
	if flag;
		flag = flag && logical( ~isempty(C.xf) );
	end
end %has_scan()

function [x,c] = push_camera_state(C,x,c)
	x=[x;C.t;C.a];
	c=[c C.edge(2)];
end %push_camera_state()

function [x,t] = push_scan_state(C,x,t)
	% FIXME: scan merging, not implemented yet
	x=[x;C.xf];
	t=t;
end %push_scan_state()

function [p1,p2]=get_correspondence(C,kpts,options)
	p1=kpts{C.edge(1)}(C.matches(1,:),1:2);
	p1=calibrate_image_points(p1,options,C.edge(1))';
	p2=kpts{C.edge(2)}(C.matches(2,:),1:2);
	p2=calibrate_image_points(p2,options,C.edge(2))';
end %get_correspondence()

function cr=calibrate_image_points(cr,options,k)
	[K,kc]=get_intrinsics(options,k);
	x=remove_lens_distortion(cr(:,1:2),kc,K);
	x(:,1)=(x(:,1)-K(1,3))/K(1,1);
	x(:,2)=(x(:,2)-K(2,3))/K(2,2);
	cr(:,1:2)=x(:,1:2);
end %calibrate_image_points()
