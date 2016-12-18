function C = initialise_graph_constraints(C,kpts,Pkin,options)
% [C,kpts]=build_camera_graph(options)
%
% Initialises edges in the camera graph with motion parameters 
% if (there is an edge)
%	initialise from kinematics
%	if (there is calibration)
%		initialise from calibration data
%	elseif (there is enough correspondence)
%		initialise from vision
%
% Tariq Abuhashim - 2016
%
% iCub - Koroibot

config_visual;

for k=1:length(C)
	% get the poeses
	if is_calibrated(C(k),options)
		C(k).t=options.t;
		C(k).a=options.a;
	elseif is_correspondence(C(k),options)&&USE_VISION
		% currently, its better not to use the essential matrix
	else
		ci=camera_matrix_to_pose(Pkin{C(k).edge(1)});
		cj=camera_matrix_to_pose(Pkin{C(k).edge(2)});
		cji=transform_to_relative_w(cj,ci);
		C(k).t=cji(1:3);
		C(k).a=cji(4:6);
	end
	% get the scans
	if is_scan(C(k),options)
		p1=kpts{C(k).edge(1)}(C(k).matches(1,:),1:2);
		p1=calibrate_image_points(p1,options,k)';
		p2=kpts{C(k).edge(2)}(C(k).matches(2,:),1:2);
		p2=calibrate_image_points(p2,options,k)';
		xs=[C(k).t;C(k).a];
		C(k).xf=test_triangulate_inverse_depth(p1,p2,xs);
	end
end
end %initialise_graph_constraints()
%
%
function flag=is_calibrated(C,options)
flag=logical( mod(C.edge(1),2) & (C.edge(2)-C.edge(1))==1 & isfield(options,'t') & isfield(options,'a') ) ;
end %is_calibrated()

function flag=is_correspondence(C,options)
flag=logical( C.weight>options.mincorrnr: ) ;
end %is_correspondence()

function flag=is_scan(C,options)
flag=logical( mod(C.edge(1),2) & is_correspondence(C,options) & (C.edge(2)-C.edge(1))==1 ) ;
end %is_scan()

function x=camera_matrix_to_pose(P)
x=[P(1:3,4);R2w(P(1:3,1:3))'];
end %camera_matrix_to_pose()

function cr=calibrate_image_points(cr,options,k)
[K,kc]=get_intrinsics(options,k); % remove lense distortion from the key-points
x=remove_lens_distortion(cr(:,1:2),kc,K);
x(:,1)=(x(:,1)-K(1,3))/K(1,1);
x(:,2)=(x(:,2)-K(2,3))/K(2,2);
cr(:,1:2)=x(:,1:2);
end %calibrate_image_points()
