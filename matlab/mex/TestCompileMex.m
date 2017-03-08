nrho = 500; % number of inverse depth data
ncams = 20; % number of cameras
C = [];
for i = 1:ncams*nrho;
    C(i).cam = randi(ncams);
    C(i).kpt = randi(nrho);
    C(i).p1 = rand(2,1); % p1
    C(i).z = rand(2,1); % p2
    C(i).R = 20*eye(2);
end
N = ncams*6 + nrho;
x = abs(rand(N,1));
xs = abs(rand(N,1));
addpath('../batch/estimator_relative');
addpath('../common/');
addpath('../batch/common/');
fprintf('matlab : ');
 tic;
C1 = generate_image_constraints_info_inverse_depth(C, xs, ncams);
toc
fprintf('C++ : ');
tic;
C2 = mex_generate_constraints_info_Mviews(C, xs, ncams);
toc