% a consule to compile and test mex files under construction
clc; clear;
cd('~/Dropbox/code/matlab/cpp/mex/');
addpath('../../batch/common/');
addpath('../../common/');
addpath('../../sba/code');
run('../../icub/set_folders');
warning off MATLAB:mex:GccVersion_link

if 1
    % Bundle-Adjustment case:
    addpath('../../batch/estimator_relative_inverse_depth_Mviews/');
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
    P = abs(rand(N));
    P = (P+P')/2;
    
    %mex mex_constraints_addition_inverse_depth_Mviews.cpp -largeArrayDims;
    suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
    mex('-v', '-largeArrayDims','mex_optimise_constraints_Mviews.cpp',...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lmwlapack','-lmwblas','-lblas');
    
    for i=1:100
        tic;
        mex_optimise_constraints_Mviews(C, xs, ncams);
        toc
    end
    
    return
    
    %     fprintf('matlab : ');
    %     tic;
    %     C1 = generate_image_constraints_info_inverse_depth(C, xs, ncams);
    %     toc
    %     mex mex_generate_constraints_info_Mviews.cpp -largeArrayDims;
    %     fprintf('C++ : ');
    %     tic;
    %     C2 = mex_generate_constraints_info_Mviews(C, xs, ncams);
    %     toc
    
    %     fprintf('matlab : ');
    %     tic;
    %     N = length(xs)/6*ncams;
    %     [yon1, Yon1] = update_info_matrix_inverse_depth(C2, N, ncams);
    %     toc
    %     mex mex_update_info_matrix_Mviews.cpp -largeArrayDims;
    %     fprintf('C++ : ');
    %     tic;
    %     N = length(xs)-6*ncams;
    %     [yon2, Yon2] = mex_update_info_matrix_Mviews(C2, N, ncams);
    %     toc
    %     full(yon1-yon2)
    %     full(Yon1-Yon2)
    
    %         fprintf('matlab : ');
    %         tic;
    %         gate_1 = compute_gate_inverse_depth(x, P, C, xs, ncams);
    %         toc
    %         mex mex_compute_gate_inverse_depth_Mviews.cpp -largeArrayDims;
    %         fprintf('C++ : ');
    %         tic;
    %         gate_2 = mex_compute_gate_inverse_depth_Mviews(x, P, C, xs, ncams);
    %         toc
    %         % number_of_gate_errors = sum(abs(gate_2-gate_1)>1e-3)
    %         error = abs((gate_2-gate_1)./gate_1) > 1e-5;
    %         number_of_gate_errors = sum(error)
    
    %     Y = 100*speye(length(xs));
    %     y = rand(length(xs), 1);
    %     C2 = mex_generate_constraints_info_Mviews(C, xs, ncams);
    %     sw = zeros(1,length(C2));
    %     options.gateinnov = chi_square_bound(.99, 2);
    %     options.verbose = 1;
    %     fprintf('matlab : ');
    %     tic;
    %     [y1, Y1, sw] = constraints_addition_inverse_depth(y, Y, C2, sw, xs, options, ncams);
    %     toc
    %     %mex mex_constraints_addition_inverse_depth_Mviews.cpp -largeArrayDims;
    %     suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
    %     mex('-v', '-largeArrayDims',...
    %         ['-I' suitesparse '/UFconfig'],...
    %         ['-I' suitesparse '/CHOLMOD/Include'],...
    %         ['-I' suitesparse '/CHOLMOD/MATLAB'],...
    %         ['-I' suitesparse '/AMD/Include'],...
    %         ['-I' suitesparse '/COLAMD/Include'],...
    %         ['-I' suitesparse '/CCOLAMD/Include'],...
    %         ['-I' suitesparse '/CAMD/Include'],...
    %         ['-I' suitesparse '/metis-4.0/Lib'],...
    %         'mex_constraints_addition_inverse_depth_Mviews.cpp',...
    %         '-lcholmod','-lmwlapack','-lmwblas','-lblas');
    %
    %     fprintf('C++ : ');
    %     tic;
    %     [y2, Y2] = mex_constraints_addition_inverse_depth_Mviews(y, Y, C2(sw==1), xs, ncams);
    %     toc
    %     %y1-y2
    %     %full(Y1-Y2)
    
    
    %     load airfoil;
    %     % Scaling x and y
    %     x = pow2(x,-32);
    %     y = pow2(y,-32);
    %     % Forming the sparse adjacency matrix and making it positive definite
    %     n = max(max(i),max(j));
    %     A = sparse(i,j,-1,n,n);
    %     A = A + A';
    %     d = abs(sum(A)) + 1;
    %     A = A + diag(sparse(d));
    %     m = symamd(A); % 'Approximate minimum degree'
    %     Y = A(m, m);
    %     y = rand(n,1);
    %
    %     %Y = speye(10);
    %     %Y(1,5)=1; Y(3,6)=1; Y(4,5)=1;
    %     %Y = (Y+Y')/2;
    %     %y = rand(10,1);
    %
    %     %fprintf('matlab : ');
    %     %tic;
    %     %[x1, P1] = recover_moments(b, A);
    %     %toc
    %     suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
    %     mex('-v', '-largeArrayDims',...
    %         ['-I' suitesparse '/UFconfig'],...
    %         ['-I' suitesparse '/CHOLMOD/Include'],...
    %         ['-I' suitesparse '/CHOLMOD/MATLAB'],...
    %         ['-I' suitesparse '/AMD/Include'],...
    %         ['-I' suitesparse '/COLAMD/Include'],...
    %         ['-I' suitesparse '/CCOLAMD/Include'],...
    %         ['-I' suitesparse '/CAMD/Include'],...
    %         ['-I' suitesparse '/metis-4.0/Lib'],...
    %         'mex_recover_moments.cpp',...
    %         '-lcholmod','-lmwlapack','-lmwblas','-lblas');
    %     %fprintf('C++ : ');
    %
    %     tic;
    %     [x2, P2] = mex_recover_moments(Y, y);
    %     toc
    %
    %     tic;
    %         P3 = spinv(Y);
    %         i = amd(Y); % compute reordering
    %         L = lchol(Y(i,i));
    %         f = L\y(i,1); % triangular solve: L\y
    %         x3(i,1) = full(L'\f); % mean with original ordering
    %         toc
    
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

% suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
% file = 'mex_recover_moments.cpp';
% d = '-largeArrayDims ';
% %d = [d ' -DLONG -D''LONGBLAS=UF_long'' -DDLONG'];
% include = ['-I' suitesparse '/AMD/Include ' ...
%     '-I' suitesparse '/CAMD/Include '...
%     '-I' suitesparse '/CCOLAMD/Include '...
%     '-I' suitesparse '/CHOLMOD/Include '...
%     '-I' suitesparse '/COLAMD/Include ' ...'
%     '-I' suitesparse '/UFconfig'];
% include = [include ' -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE'];
% metis_path = '../../metis-4.0' ;
% include = [include ' -I' metis_path '/Lib'];
% extra_path = '/home/tabuhashim/Dev/cholmod-extra';
% include = [include ' -I' extra_path '/Include'];
% lapack = ' ';
% %lapack = [lapack '-lmwlapack -lmwblas'];
% source = ['-I' suitesparse '/CHOLMOD/Lib'];
% s = sprintf('mex -c %s %s %s %s %s', file, include, source, d, lapack);
% fprintf('%s\n',s);
% eval(s);

% mex mex_recover_moments.cpp ...
%     -largeArrayDims ...
%     -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/Include ...
%     -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/CHOLMOD/MATLAB ...
%     -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/UFconfig ...
%     -I/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse/metis-4.0/Lib ...
%     -I/home/tabuhashim/Dev/cholmod-extra/Include;

% suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
% mex('mex_recover_moments.cpp',...
%     ['-I' suitesparse '/UFconfig'],...
%     ['-I' suitesparse '/CHOLMOD/Include'],...
%     ['-I' suitesparse '/CHOLMOD/MATLAB'],...
%     ['-I' suitesparse '/metis-4.0/Lib'],...
%     '-lcholmod','-lcholmod','-largeArrayDims');