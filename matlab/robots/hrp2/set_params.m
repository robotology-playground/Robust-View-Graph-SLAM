function options = set_params(options, field, val)
%options = set_params(options)
% loads options for the batch implementation
%
%	For HRP2@Toulouse.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Auckland, New Zealand, 2017

if nargin < 1; % set basic defaults
    options.first_image = 0;
    options.last_image = 0;
    options.steps = 0;
    options.verbose = 2;
    return
end

if nargin == 1 % set default vision and optimiser configurations
    options = set_vision_options(options);
    options = set_optmsr_options(options);
end

if nargin == 2 % assign calibration or vision or optimiser configurations
    if strcmp(field, 'calib')
        options = set_calib_options(options);
    elseif strcmp(field, 'vision')
        options = set_vision_options(options);
    elseif strcmp(field, 'optimiser')
        options = set_optmsr_options(options);
    else
        error('Invalid option, should be "vision" or "optimiser"')
    end
end

if nargin == 3 % change a field
    if strcmp(field, 'folder')
        options.(field) = val;
    elseif strcmp(field, 'save')
        options.(field) = val;
    elseif strcmp(field, 'calib')
        options.(field) = val;
        options = set_calib_options(options);
    else
        assert(isfield(options, field), ['Invalid field ' field]);
        options.(field) = val;
    end
end

if isfield(options, 'save')
    if ~isdir(strcat(options.save))
        mkdir(strcat(options.save));
    end
    save(strcat(options.save,'/options'), 'options');
end

function options = set_vision_options(options)
% ============= % vision options - defaults
options.splitimage  =0;
%options.minbase     =34/1000;% minimum baseline to perform triangulation (meters)
options.mindisp     =2;% minimum pixel displacement in the images to be considered
%options.bucketsize      =[150,150];% bucket size
%options.ransac          =200;% number of ransac iterations
%options.RANSAC_pixtol   =.5;% tolerance RANSAC   .1 with ladybug / 1 marulan
options.mincorrnr       =25;% minimum number of inliers to trust two-view results
options.mininlnr        =10;% min num of matches to compute two-view geometries
%options.merge_tracks=1;% perform global data associations?
%options.roterrtol=eps;% threshold in pixels for inlier (in rotation averaging)
%options.maxnumfeat      =Inf;
%options.gridsize        =0;% 3
options.gridmargin      =5;% 20
options.gridhorizon     =5;% 120
options.detector  ='KAZE';% sift, kaze
if strcmp(options.detector,'SIFT')
    options.siftthreshold   =[50,0];
    options.siftratio       =2;
elseif strcmp(options.detector,'KAZE')
    kazethreshold   =0.0001;
    kazeratio       =0.5;
    options.detector_param =[kazethreshold,kazeratio];
elseif strcmp(options.detector,'FAST')
    fastmargin      =10;%20
    fastthreshold   =20;% 20
    fastnonmax      =0;% 0
    LKWinSize       =[15,15];
    options.detector_param =[fastmargin,fastthreshold,fastnonmax,LKWinSize(1)];
else
    error('Features extraction method can be either sift, kaze or fast');
end
% % ============= % old rotation averaging options
% opt.visviews = 10; % points should be seen in at least this many cameras
%             % (in patch_solve_known_rotation_1 and patch_solve_known_rotation_2)
% opt.pixtol1 = 2; % tolerans for first known-rotation run
%             % (in patch_solve_known_rotation_1)
% opt.pixtol2 = 2; % tolerans for second known-rotation run
%             % (in patch_solve_known_rotation_2)
% opt.middepth = .1; % minimum depth bound
%             % (in patch_solve_known_rotation_1 and patch_solve_known_rotation_2)
% opt.maxdepth = 5; % maximum depth bound
%             % (in patch_solve_known_rotation_1 and patch_solve_known_rotation_2)
% opt.uncertin_tol = 10; % depth uncertainty to accept points in reconstruction
%
%
function options = set_optmsr_options(options)  
% ============= % epipolar constraints optimiser options - defaults
options.ncams       =20;% number of frames inside the bundle
options.nkeys       =20;% number of key-frames in the bundle
options.nview       =15;% minimum number of views to consider a point in the map
%options.sigma_a    =1;% rotational noise std in degrees
%options.sigma_b    =1;% translational noise std in milli-meters
%options.sigma_c    =.5;% pixels (3d points uncertainty)
options.sigma_r     =0.5;% image noise std in pixels
options.gateratio   =0.5;% 0.5 - larger is slower
options.gateinnov=chi_square_bound(.99,2);% innovations gate (default)
options.gateresid=chi_square_bound(.95,2);% residuals gate (default)
options.gatetrust=chi_square_bound(.95,2);% trusted constraints gate (default)
options.iteration   =2;
options.checkrank   =1;
options.verbose     =2;% verbose level: 0,1,2
options.maxitr      =100;% number of iterations, and residuals gate thresholds
%
%
function options = set_calib_options(options)
% ============= % Camera calibration 
out=read_calibration_file(options.calib);
options.imgsize=out([1 2]);
options.K1=[out(3),0,out(5);0,out(4),out(6);0,0,1];
options.fc1=out([3 4]);
options.cc1=out([5 6]);
options.kc1=out([23:26]);
options.K2=[out(7),0,out(9);0,out(8),out(10);0,0,1];
options.fc2=out([7 8]);
options.cc2=out([9 10]);
options.kc2=out([27:30]);
options.R=reshape(out([11:19]),3,3);
options.t=out([20:22]);
