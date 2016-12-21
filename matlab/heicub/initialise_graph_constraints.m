function C = initialise_graph_constraints(C,kpts,Pkin,options)
%[C,kpts]=build_camera_graph(options)
%
% Initialises edges in the camera graph with motion parameters 
% if (there is an edge)
%	if (there is calibration)
%		initialise from calibration data
%	elseif (there is enough correspondence)
%		initialise from vision
%	elseif (there are kinematics)
%		initialise from kinematics
%	else
%		random initialisation
%
% Tariq Abuhashim - 2016
% t.abuhashim@gmail.com
%
% iCub - Koroibot

config_visual;

for k=1:length(C)
	% get the poeses
	if is_calibrated(C(k),options)
		C(k).t=options.t;
		C(k).a=options.a;
	elseif has_correspondence(C(k),options)&&USE_VISION
		% currently, its better not to use the essential matrix
	elseif has_kinematics(Pkin{C(k).edge})
		ci=camera_matrix_to_pose(Pkin{C(k).edge(1)});
		cj=camera_matrix_to_pose(Pkin{C(k).edge(2)});
		cji=transform_to_relative_w(cj,ci);
		C(k).t=cji(1:3);
		C(k).a=cji(4:6);
	else
		C(k).t=rand(3,1);
		C(k).a=rand(3,1);
		warning(['Edge ', num2str(C(k).edge(1)),'-',num2str(C(k).edge(2)), ... 
				' could not be initialised and is set to RAND']);	
	end
	% get the scans
	if has_scan(C(k),options)
		[p1,p2]=get_correspondence(C(k),kpts,options);
		xs=[C(k).t;C(k).a];
		C(k).xf=test_triangulate_inverse_depth(p1,p2,xs)';
	end
end
end %initialise_graph_constraints()

function flag=is_calibrated(C,options)
	flag=logical( mod(C.edge(1),2) & (C.edge(2)-C.edge(1))==1 ...
		& isfield(options,'t') & isfield(options,'a') ) ;
end %is_calibrated()

function flag=has_correspondence(C,options)
	flag=logical( C.weight>options.mincorrnr ) ;
end %is_correspondence()

function flag=has_kinematics(P1,P2)
	flag=logical( ~isempty(P1) & ~isempty(P2) );
end

function flag=has_scan(C,options)
	% FIXME: scan merging, not implemented yet
	flag=logical( mod(C.edge(1),2) & has_correspondence(C,options) ... 
		& (C.edge(2)-C.edge(1))==1 ) ;
end %is_scan()

function x=camera_matrix_to_pose(P)
	x=[P(1:3,4);R2w(P(1:3,1:3))'];
end %camera_matrix_to_pose()

function [p1,p2]=get_correspondence(C,kpts,options)
	p1=kpts{C.edge(1)}(C.matches(1,:),1:2);
	p1=calibrate_image_points(p1,options,C.edge(1))';
	p2=kpts{C.edge(2)}(C.matches(2,:),1:2);
	p2=calibrate_image_points(p2,options,C.edge(2))';
end %get_correspondence()

function cr=calibrate_image_points(cr,options,k)
	[K,kc]=get_intrinsics(options,k); % remove lense distortion from the key-points
	x=remove_lens_distortion(cr(:,1:2),kc,K);
	x(:,1)=(x(:,1)-K(1,3))/K(1,1);
	x(:,2)=(x(:,2)-K(2,3))/K(2,2);
	cr(:,1:2)=x(:,1:2);
end %calibrate_image_points()
