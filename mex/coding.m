% a consule to compile and test mex files under construction
clc; clear;
cd('~/Dropbox/code/matlab/cpp/mex/');
addpath('../../batch/common/');
addpath('../../common/');
addpath('../../sba/code');
warning off MATLAB:mex:GccVersion_link

if 1
    % Bundle-Adjustment case:
    addpath('../../batch/estimator_relative_inverse_depth_Mviews/');
    
    nrho = 1; % number of inverse depth data
    ncams = 2; % number of cameras
    C = [];
    for i = 1:ncams*nrho;
        C(i).cam = randi(ncams);
        C(i).kpt = randi(nrho);
        C(i).p1 = rand(2,1); % p1
        C(i).z = rand(2,1); % p2
        C(i).R = eye(2);
    end
    N = ncams*6 + nrho;
    x = abs(rand(N,1));
    xs = abs(rand(N,1));
    P = abs(rand(N));
    P = (P+P')/2;
    
    fprintf('matlab : ');
    tic;
    C1 = generate_image_constraints_info_inverse_depth(C, xs, ncams);
    toc
    mex mex_generate_constraints_info_Mviews.cpp -largeArrayDims;
    fprintf('C++ : ');
    tic;
    C2 = mex_generate_constraints_info_Mviews(C, xs, ncams);
    toc
    
%     fprintf('matlab : ');
%     tic;
%     [yon1, Yon1] = update_info_matrix_inverse_depth(C1, N, ncams);
%     toc
%     warning('the function mex_update_info_matrix_Mviews is not implemented yet'); 
%     mex mex_update_info_matrix_Mviews.cpp -largeArrayDims;
%     fprintf('C++ : ');
%     tic;
%     [yon2, Yon2] = mex_update_info_matrix_Mviews(C2, N, ncams);
%     toc
    
    %     fprintf('matlab : ');
    %     tic;
    %     gate_1 = compute_gate_inverse_depth(x, P, C, xs, ncams);
    %     toc
    %     mex mex_compute_gate_inverse_depth_Mviews.cpp -largeArrayDims;
    %     fprintf('C++ : ');
    %     tic;
    %     gate_2 = mex_compute_gate_inverse_depth_Mviews(x, P, C, xs, ncams);
    %     toc
    %     % number_of_gate_errors = sum(abs(gate_2-gate_1)>1e-3)
    %     error = abs((gate_2-gate_1)./gate_1) > 1e-5;
    %     number_of_gate_errors = sum(error)
    
else
    % Pose-graph case:
    addpath('../../batch/estimator_global/');
    mex mex_compute_gate_graph.cpp -largeArrayDims;
    ncams = 100;
    ncon = ncams*(ncams-1)/2; % dense graph
    for i = 1:ncon;
        C(i).edge = [randi(ncams) randi(ncams)];
        %C(i).edge = [1 2];
        C(i).z = rand(6,1);
        C(i).R = 2*eye(6);
    end
    x = rand(ncams*6,1);
    xs = rand(ncams*6,1);
    P = abs(rand(ncams*6));
    P = (P+P')/2;
    fprintf('matlab : ');
    tic;
    gate_1 = compute_residuals(x, P, C, xs);
    toc
    fprintf('C++ : ');
    tic;
    gate_2 = mex_compute_gate_graph(x, P, C, xs);
    toc
    
    error = abs((gate_2-gate_1)./gate_1) > 1e-2;
    number_of_constraints = ncon
    number_of_gate_errors = sum(error)
    gate_1(error==1)
    gate_2(error==1)
end