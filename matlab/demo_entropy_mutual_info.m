% demo on Entropy and Mutual Information

folder = '/home/tabuhashim/Documents/data/WalkingDatasets/results/run_2/';
run('../icub/set_folders');
load(strcat(folder,'options'));
load(strcat(folder,'encoders'));
load(strcat(folder,'floatingbase'));
load(strcat(folder,'constraints'));

k = 1; % reference image
ncams = 20; % number of frames in the bundle

% get kinematics
[Pkin, A0] = cameras_from_kinematics(eyes, neck, waist, floatingbase);

% get point tracks
p = get_aligned_point_matches_v2(options, k, ncams);
ncams = length(p) ;% final number of images (in case some were removed)
npnts = length(p(1).x) ;% number of points tracked, the same for all images

% initialise image constraints
[C, Ct, sw] = initialise_constraints(p, options);

% initialistion point
xs = initialise_linearisation_point_v2(options, p, Pkin);

% information initialisation
if ~isempty(Ct)
    Ct = mex_generate_constraints_info_Mviews(Ct, xs, ncams);
end
[y, Y] = initialise_info_matrix(Ct, xs, ncams); % sw is for C, not for Ct

% generate constraints information
C = mex_generate_constraints_info_Mviews(C, xs, ncams);

% compute information measures for the second frame
k = 8;
H = zeros(npnts,1); % Entropy
M = zeros(npnts,1); % Mutual information
for i = 1 : length(C)
    if C(i).cam == k && any(C(i).Y(:))
        H(C(i).kpt,1) = trace(C(i).Y);
        I(C(i).kpt,1) = 0;
    end
end