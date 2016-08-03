% a consule to compile and test mex files under construction
clc; clear;
cd('~/Dropbox/code/matlab/cpp/mex/');
%addpath('../../batch/common/');
%addpath('../../common/');
%addpath('../../sba/code');
run('../../icub/set_folders');
warning off MATLAB:mex:GccVersion_link

load airfoil;
% Scaling x and y
x = pow2(x,-32);
y = pow2(y,-32);
% Forming the sparse adjacency matrix and making it positive definite
n = max(max(i),max(j));
A = sparse(i,j,-1,n,n);
A = A + A';
d = abs(sum(A)) + 1;
A = A + diag(sparse(d));
m = symamd(A); % 'Approximate minimum degree'
Y = A(m, m);
y = rand(n,1);

%Y = speye(10);
%Y(1,5)=1; Y(3,6)=1; Y(4,5)=1;
%Y = (Y+Y')/2;
%y = rand(10,1);

%fprintf('matlab : ');
%tic;
%[x1, P1] = recover_moments(b, A);
%toc

suitesparse = '/home/tabuhashim/Dev/GPstuff-4.6/SuiteSparse';
mex('-v', '-largeArrayDims',...
    'mex_recover_moments.cpp',...
    ['-I' suitesparse '/UFconfig'], ['-I' suitesparse '/CHOLMOD/Include'],...
    ['-I' suitesparse '/CHOLMOD/MATLAB'], ['-I' suitesparse '/AMD/Include'],...
    ['-I' suitesparse '/COLAMD/Include'], ['-I' suitesparse '/CCOLAMD/Include'],...
    ['-I' suitesparse '/CAMD/Include'], ['-I' suitesparse '/metis-4.0/Lib'],...
    '../src/RecoverMoments.cpp',...
    '-lcholmod','-lmwlapack','-lmwblas','-lblas') ;

tic;
P1 = spinv(Y);
i = amd(Y); % compute reordering
L = lchol(Y(i,i));
f = L\y(i,1); % triangular solve: L\y
x1(i,1) = full(L'\f); % mean with original ordering
fprintf('Matlab: ');toc

tic;
[x2, P2] = mex_recover_moments(Y, y);
fprintf('C++: ');toc

spy(P2-P1); title('patters on P2-P1');
fprintf(['Solving error = ', num2str(sum(norm(x2-x1))),'\n']);
fprintf(['Inverse error = ', num2str(sum(norm(P2(:)-P1(:)))),'\n']);