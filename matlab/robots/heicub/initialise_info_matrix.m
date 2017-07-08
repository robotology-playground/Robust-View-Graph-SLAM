function [y,Y]=initialise_info_matrix(Ct,xs,ncams,options)
%[y,Y]=initialise_info_matrix(Ct,xs,ncams)
%
%	Information matrix initialisation for iCub@heidelberg
%	This assumes the calibrated stereo case.
%
%	For iCub@heidelberg.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014

config_visual();

s0 = 0.00001;
if is_calibrated(xs(7:12),options)
    t2 = options.t_error/1000; % meters, c=2, from options.t_error
    a2 = options.a_error*180/pi; % degrees, c=2, from options.r_error
end
if USE_VISION
    t3 = 0.0005; % meters, c=3:2:end
    a3 = 1; % degrees, c=3:2:end
    t4 = 0.0005; % meters, c=4:2:end
    a4 = 1; % degrees, c=4:2:end
    sd = 1; % depth (meters)
else
    t3 = 10/1000; % meters, c=3:2:end
    a3 = 2; % degrees, c=3:2:end
    t4 = 10/1000; % meters, c=4:2:end
    a4 = 2; % degrees, c=4:2:end
    sd = 0.5; % depth (meters)
end

% consider trusted constraints in Ct, if empty, initialise empty {y,Y}
if ~isempty(Ct) % Generate trusted information
    Ct = mex_generate_constraints_info_Mviews(Ct, xs, ncams);
end
N = length(xs);
npts = N-6*ncams;
[y, Y] = update_info_matrix_inverse_depth(Ct, npts, ncams);

% prior information (Easy initialisation, no uncertainty modeling yet)
s = ones(N,1)./(sd^2); % depth error ( .5 meters with vision->R , 1 meters with kinematics->R )
s(1:6) = 1/(s0)^2; % cam 1, 0.0001 error in global
for j = 1:ncams
	if is_calibrated(xs(7:12),options)
        s((j-1)*6+(1:3)) = 1./(t2).^2; % .0005 meters
        s((j-1)*6+(4:6)) = 1./(a2.*pi/180).^2; % .5 degrees
        continue
    end
    if mod(j,2)==1
		s((j-1)*6+(1:3)) = 1/(t3)^2; % .0005 meters
    	s((j-1)*6+(4:6)) = 1/(a3*pi/180)^2; % .5 degrees with vision
    else
    	s((j-1)*6+(1:3)) = 1/(t4)^2; % .0005 meters
    	s((j-1)*6+(4:6)) = 1/(a4*pi/180)^2; % .5 degrees with vision
    end
end

% update information matrix and vector
Y0 = sparse(1:N, 1:N, s, N, N);
y = y + Y0*xs;
Y = Y + Y0;

% Y = speye(N, N);
% y = .01*rand(N, 1);

end

% Additional functions
function flag=is_calibrated(x,options)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
	flag=logical( isfield(options,'t') & isfield(options,'a') );
	if flag
		t=-w2R(options.a)*options.t/1000; % also in initialise_graph_constraints.m
		a=-options.a;
		flag=logical( x(1:3)==t & x(4:6)==a ) ;
	end
end %is_calibrated()
