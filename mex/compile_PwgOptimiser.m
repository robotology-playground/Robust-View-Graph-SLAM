% a consule to compile and test mex files under construction
clc; clear;
warning off MATLAB:mex:GccVersion_link

vgslam = '/home/tariq/Documents/Robust-View-Graph-SLAM';
eigen = '/home/tariq/Downloads/Eigen' ;
suitesparse = '/home/tariq/Downloads/SuiteSparse' ;

cd([vgslam,'/mex/']);
run('../../heicub/set_folders');

model = 1 ;
jacobian = 1 ;
generate = 1 ;
update = 1 ; if (update&&~generate); generate = 1 ; end;
gate = 1 ;
addition = 0 ; if (addition&&~generate); generate = 1 ; end;
optimise = 0 ; if (optimise&&~generate); generate = 1 ; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (model)
    mex('-v', '-largeArrayDims',...
        'mex_observation_model_inverse_depth.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
    
    mex('-v', '-largeArrayDims',...
        'mex_observation_model_jacobian_inverse_depth_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (jacobian)
    mex('-v', '-largeArrayDims',...
        'mex_observation_model_jacobian_inverse_depth.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
    
    mex('-v', '-largeArrayDims',...
        'mex_observation_model_inverse_depth_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (generate)
    mex('-v', '-largeArrayDims',...
        'mex_generate_constraints_info_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (update)
    mex('-v', '-largeArrayDims',...
        'mex_update_info_matrix_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (gate)
    mex('-v', '-largeArrayDims',...
        'mex_compute_gate_inverse_depth_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (addition)
    mex('-v', '-largeArrayDims',...
        'mex_constraints_addition_inverse_depth_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

if (optimise)
    mex('-v', '-largeArrayDims',...
        'mex_optimise_constraints_Mviews.cpp',...
        ['-I' vgslam '/include'], ...
        ['-I' eigen ], ...
        ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
        ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
        ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
        ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
        '../src/PwgOptimiser.cpp', '../src/RecoverMoments.cpp',...
        '-lcholmod','-lcolamd','-lcamd','-lamd','-lccolamd','-lmwlapack','-lmwblas','-lblas','-lmetis',['-L' suitesparse '/metis-4.0']) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%

nrho = 3000 ; % number of inverse depth data
ncams = 40 ; % number of cameras
C = [] ; k = 1 ;
for i = 1:ncams
    for j = 1:nrho
        C(k).cam = i ;
        C(k).kpt = j ;
        C(k).p1 = rand(2,1) ; % p1
        C(k).z = C(k).p1+.005 ; % p2
        C(k).R = eye(2) ;
        k = k + 1 ;
    end
end
N = ncams*6 + nrho ;
x(1:6,1) = 0;
for i = 2:ncams
    x((i-1)*6+(1:3),1) = [.068; 0; 0];
    x((i-1)*6+(1:3),1) = [0; 0; 0];
end
x(6*ncams+(1:nrho),1) = .2*rand(nrho,1) + .5;
xs = x + .1*rand(6*ncams+nrho,1) + .01;

%%%

if (generate)
    fprintf('\n') ;
    fprintf('####### testing generation of constraints info ..... \n');
    fprintf('\n') ;
    tic;
    C1 = generate_image_constraints_info_inverse_depth(C, xs, ncams);
    fprintf('Matlab: '); toc
    tic;
    C2 = mex_generate_constraints_info_Mviews(C, xs, ncams);
    fprintf('C++: '); toc
    cam_error = 0 ; kpt_error = 0 ; p1_error = 0 ; z_error = 0 ;
    R_error = 0 ; y_error = 0 ; Y_error = 0 ;
    eps1 = 1e-6; eps2 = 1e-3;
    for k = 1 : length(C1)
        cam_error = (C1(k).cam ~= C2(k).cam) + cam_error ;
        cam_error = (C1(k).kpt ~= C2(k).kpt) + cam_error ;
        p1_error = (max(abs(C1(k).p1  - C2(k).p1  ))>eps1) + p1_error ;
        z_error =  (max(abs(C1(k).z   - C2(k).z   ))>eps1) + z_error ;
        R_error =  (max(abs(C1(k).R(:)- C2(k).R(:)))>eps1) + R_error ;
        y_error =  (max(abs(C1(k).y   - C2(k).y   ))>eps2) + y_error ;
        Y_error =  (max(abs(C1(k).Y(:)- C2(k).Y(:)))>eps2) + Y_error ;
    end
    fprintf(['number_of_C.cam_errors = ',num2str(sum(cam_error)),'\n']);
    fprintf(['number_of_C.kpt_errors = ',num2str(sum(kpt_error)),'\n']);
    fprintf(['number_of_C.p1_errors = ',num2str(sum(p1_error)),' (eps=',num2str(eps1),')','\n']);
    fprintf(['number_of_C.z_errors = ',num2str(sum(z_error)),' (eps=',num2str(eps1),')','\n']);
    fprintf(['number_of_C.R_errors = ',num2str(sum(R_error)),' (eps=',num2str(eps1),')','\n']);
    fprintf(['number_of_C.y_errors = ',num2str(sum(y_error)),' (eps=',num2str(eps2),')','\n']);
    fprintf(['number_of_C.Y_errors = ',num2str(sum(Y_error)),' (eps=',num2str(eps2),')','\n']);
    fprintf('\n');
end

%%%

P = abs(rand(N)) ;
P = (P+P')/2 ;

%%%

if (gate)
    fprintf('\n');
    fprintf('####### testing gate computations ..... \n');
    fprintf('\n') ;
    tic;
    gate_1 = compute_gate_inverse_depth(x, P, C, xs, ncams);
    fprintf('Matlab: '); toc
    tic;
    gate_2 = mex_compute_gate_inverse_depth_Mviews(x, P, C, xs, ncams);
    fprintf('C++: '); toc
    %number_of_gate_errors = sum(abs(gate_2-gate_1)>1e-3)
    eps1 = 1e-6;
    error = abs((gate_2-gate_1)./gate_1) > eps1;
    fprintf(['number_of_gate_errors = ',num2str(sum(error)),' (eps=',num2str(eps1),')','\n']);
    fprintf('\n');
end

%%%

Y = speye(length(xs));
%Y = sparse(1:length(xs),1:length(xs),rand(length(xs),1),length(xs),length(xs));
y = zeros(length(xs),1);
sw = zeros(1,length(C2));
options.gateinnov = chi_square_bound(.99,2);
options.verbose = 1;

%%%

if (update)
    fprintf('\n') ;
    fprintf('####### testing information update ..... \n') ;
    fprintf('\n') ;
    tic ;
    N = length(xs)-(6*ncams) ;
    [yon1, Yon1] = update_info_matrix_inverse_depth(C2, N, ncams) ;
    fprintf('Matlab: ') ; toc
    tic ;
    [yon2, Yon2] = mex_update_info_matrix_Mviews(C2, N, ncams) ;
    fprintf('C++: ') ; toc
    a = abs(full(yon1-yon2)) ;
    A = abs(full(Yon1-Yon2)) ;
    eps1 = 1e-6;
    fprintf(['number_of_yon_errors = ',num2str(sum(a(:)>eps1)),' (eps=',num2str(eps1),')','\n']) ;
    fprintf(['max_yon_error = ',num2str(max(a(:))),'\n']) ;
    fprintf(['number_of_Yon_errors = ',num2str(sum(A(:)>eps1)),' (eps=',num2str(eps1),')','\n']) ;
    fprintf(['max_Yon_error = ',num2str(max(A(:))),'\n']) ;
    fprintf('\n') ;
end

%%%

if (addition)
    fprintf('\n');
    fprintf('####### testing constraints addition ..... \n');
    fprintf('\n') ;
    tic;
    [y1, Y1, sw1] = constraints_addition_inverse_depth(y, Y, C2, sw, xs, options, ncams) ;
    fprintf('Matlab: '); toc
    tic;
    [y2, Y2, sw2] = mex_constraints_addition_inverse_depth_Mviews(y, Y, C2, sw, xs, ncams) ;
    fprintf('C++: '); toc
    a = abs(full(y1-y2)) ;
    A = abs(full(Y1-Y2)) ;
    s = abs(full(sw1-sw2)) ;
    eps1 = 1e-6;
    fprintf(['number_of_y_errors = ',num2str(sum(a(:)>eps1)),'\n']) ;
    fprintf(['max_y_error = ',num2str(max(a(:))),' (eps=',num2str(eps1),')','\n']) ;
    fprintf(['number_of_Y_errors = ',num2str(sum(A(:)>eps1)),'\n']) ;
    fprintf(['max_Y_error = ',num2str(max(A(:))),' (eps=',num2str(eps1),')','\n']) ;
    fprintf(['number_of_sw_errors = ',num2str(sum(s(:)>0)),'\n']) ;
    fprintf('\n') ;
end

%%%

k = min(20,length(C)) ;
i = 6*ncams+C(k).kpt;
c = (C(k).cam-1)*6;

%%%

if (optimise)
    options.gateratio = 0.9;
    options.gateinnov = 9.2;
    options.gateresid = 9.2;
    options.gatetrust = 9.2;
    options.verbose = 0;
    Ct = [] ; sw = zeros(1, length(C2)) ;
    [y, Y, sw, x1] = optimise_constraint_image_inverse_depth( Ct, C2, sw, xs, options, ncams ) ;
    P1 = spinv(Y) ;
    [x2, P2] = mex_optimise_constraints_Mviews( C2, sw, xs, ncams ) ;
end

%%%

if (model)
    fprintf('\n') ;
    fprintf('####### testing the observation model ..... \n') ;
    fprintf('\n') ;
    tic ;
    z1 = observation_model_inverse_depth( xs([c+(1:6), i]), C(k).p1 ) ;
    fprintf('Matlab: ') ; toc
    tic ;
    z2 = mex_observation_model_inverse_depth( xs([i, c+(1:6)]), C(k).p1 ) ;
    fprintf('C++: ') ; toc
    tic ;
    z3 = mex_observation_model_inverse_depth_Mviews( xs, C(k).p1, i, c+1 ) ;
    fprintf('C++ (Mviews): ') ; toc
    fprintf(['z1 = ',num2str(z1(1)),',',num2str(z1(2)),'\n']) ;
    fprintf(['z2 = ',num2str(z2(1)),',',num2str(z2(2)),'\n']) ;
    fprintf(['z3 = ',num2str(z3(1)),',',num2str(z3(2)),'\n']) ;
    fprintf('\n') ;
end

%%%

if (jacobian)
    fprintf('\n') ;
    fprintf('####### testing the observation model jacobian ..... \n') ;
    fprintf('\n') ;
    tic ;
    H1 = observation_model_jacobian_inverse_depth( xs([c+(1:6), i]), C(k).p1 ) ;
    fprintf('Matlab: ') ; toc
    tic ;
    H2 = mex_observation_model_jacobian_inverse_depth( xs([i, c+(1:6)]), C(k).p1 ) ;
    fprintf('C++: ') ; toc
    tic ;
    H3 = mex_observation_model_jacobian_inverse_depth_Mviews( xs, C(k).p1, i, c+1 ) ;
    fprintf('C++ (Mviews): ') ; toc
    full(H1)
    full(H2(:,[2:7 1])) % inverse depth first, permutated to be consistent with H1 and H3
    full(H3)
    fprintf('\n') ;
end

%%%

fprintf('DONE !!! \n') ;
