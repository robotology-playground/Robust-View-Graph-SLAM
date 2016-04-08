clc; clear all;  close all;
run('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/icub/set_folders');
G = sparse([1 1 1 2 2 3 3 4 5 6 7 7 8 9 9  9 9], ...
    [2 6 8 3 1 4 2 5 4 7 6 4 9 8 10 5 3],1,10,10);
[s,c]=mex_SCC(full(G+G')); c=c+1; s=s+1; % cpp
m=0;j=0;for i=1:s;k=sum(c==i);if(k>m);m=k;j=i;end;end;
i=find(c==j);[u,v]=find(G+G'); I = [u v]';[~,I]=ismember(I,i);

N=max(max(I));
n=N-1;       % nodes-1 
m=size(I,2); % V

QQ=rand(m,4);
start=tic;
[Qm]=box_median(QQ,I);
%[Qm]=robust_mean(QQ,I,5,Qm);
time_m=toc(start);

start=tic;
[Qc]=mex_rotation_averaging(I-1, QQ, rand(m,3));
time_c=toc(start);

max_error=max(max([Qm-Qc]))
time_multiplier = time_m/time_c