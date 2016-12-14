function [Q,t,sptree]=demo_spanning_tree_old(G,pwg,opt)

fprintf('Performing rotation averaging ... \n');

% find the largest connected component (Graph needs to be binary)
tic;
G2=(G>0);G2=G2+G2';
[u,v]=find(G2);edges=[u v]';
[s,c,j]=mex_SCC(full(G2));
s=s+1; % cpp
c=c+1;
j=(j==1);
i=find(c==s);
[~,edges]=ismember(edges,i);  % re-number the graph after removing nodes
% which are not part of SCC
fprintf(['[RotAvg] # graph constraints = ',num2str(sum(G2(:))),'.\n']);
fprintf(['[RotAvg] # SSC constraints   = ',num2str(sum(j))    ,'.\n']);

% information weighted graph
Gy=zeros(size(G));
for ii=1:size(G,1);
    for jj=ii+1:size(G,2);
        if ~G(ii,jj);continue;end;
        %Gy(ii,jj) = length(pwg{ii,jj}.inliers);
        %Gy(ii,jj) = length(pwg{ii,jj}.maxinlier);
        Gy(ii,jj)=trace(pwg{ii,jj}.Y(4:6,4:6));
    end
end

% perform motion avaraging (replace this with robust least-squares later)
% initial guess using maximum spannig tree estimate
% pair-wise constraints
Qm=zeros(size(edges,2),4);
tm=zeros(size(edges,2),3);
for ii=1:size(edges,2)
    if edges(1,ii)<edges(2,ii)
        Qm(ii,:)=pwg{edges(1,ii),edges(2,ii)}.Q;
        tm(ii,:)=pwg{edges(1,ii),edges(2,ii)}.t;
    else
        temp=pwg{edges(2,ii),edges(1,ii)}.Q;
        R=inv(q2R(temp));
        Qm(ii,:)= R2q(R);
        tm(ii,:)=-R*pwg{edges(2,ii),edges(1,ii)}.t;
    end
end
edges=edges(:,j);
Qm=Qm(j,:);
tm=tm(j,:);

itr=0;
N=max(max(edges));
while (itr<100)
    
    [Q,t,sptree]=initialise_from_a_tree(Qm,tm,edges); % from a spannning tree
    w=zeros(N,3);
    for j=1:N
        w(j,:)=R2w(q2R(Q(j,:)));
    end
    
    % tilt
    figure(1);
    subplot(3,1,1);plot(w(1:2:end,1)*180/pi,'color',[1 .75 .75]);drawnow; axis tight;grid on;hold on;
    figure(2);
    subplot(3,1,1);plot(w(2:2:end,1)*180/pi,'color',[1 .75 .75]);drawnow; axis tight;grid on;hold on;
    % heading
    figure(1);
    subplot(3,1,2);plot(w(1:2:end,2)*180/pi,'color',[.75 1 .75]);drawnow; axis tight;grid on;hold on;
    figure(2);
    subplot(3,1,2);plot(w(2:2:end,2)*180/pi,'color',[.75 1 .75]);drawnow; axis tight;grid on;hold on;
    % roll
    figure(1);
    subplot(3,1,3);plot(w(1:2:end,3)*180/pi,'color',[.75 .75 1]);drawnow; axis tight;grid on;hold on;
    figure(2);
    subplot(3,1,3);plot(w(2:2:end,3)*180/pi,'color',[.75 .75 1]);drawnow; axis tight;grid on;hold on;
    % x
    figure(3);
    subplot(3,1,1);plot(t(1:2:end,1),'color',[1 .75 .75]);drawnow; axis tight;grid on;hold on;
    figure(4);
    subplot(3,1,1);plot(t(2:2:end,1),'color',[1 .75 .75]);drawnow; axis tight;grid on;hold on;
    % y
    figure(3);
    subplot(3,1,2);plot(t(1:2:end,2),'color',[.75 1 .75]);drawnow; axis tight;grid on;hold on;
    figure(4);
    subplot(3,1,2);plot(t(2:2:end,2),'color',[.75 1 .75]);drawnow; axis tight;grid on;hold on;
    % z
    figure(3);
    subplot(3,1,3);plot(t(1:2:end,3),'color',[.75 .75 1]);drawnow; axis tight;grid on;hold on;
    figure(4);
    subplot(3,1,3);plot(t(2:2:end,3),'color',[.75 .75 1]);drawnow; axis tight;grid on;hold on;
    
    itr=itr+1;
    
end

[Q,t,sptree]=initialise_from_MST(Qm,tm,G); % from the maximum weighted spannning tree
w=zeros(N,3);
for j=1:N
    w(j,:)=R2w(q2R(Q(j,:)));
end

% tilt
figure(1);
subplot(3,1,1);plot(w(1:2:end,1)*180/pi,'color',[1 .0 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_1$ (deg.)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(2);
subplot(3,1,1);plot(w(2:2:end,1)*180/pi,'color',[1 .0 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_1$ (deg.)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
% heading
figure(1);
subplot(3,1,2);plot(w(1:2:end,2)*180/pi,'color',[.0 1 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_2$ (deg.)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(2);
subplot(3,1,2);plot(w(2:2:end,2)*180/pi,'color',[.0 1 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_2$ (deg.)','interpreter','latex','fontsize',20);
% roll
figure(1);
subplot(3,1,3);plot(w(1:2:end,3)*180/pi,'color',[.0 .0 1],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_3$ (deg.)','interpreter','latex','fontsize',20);
xlabel('Node indices','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(2);
subplot(3,1,3);plot(w(2:2:end,3)*180/pi,'color',[.0 .0 1],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$v_3$ (deg.)','interpreter','latex','fontsize',20);
xlabel('Node indices','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
% x
figure(3);
subplot(3,1,1);plot(t(1:2:end,1),'color',[1 .0 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_1$ (mm)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(4);
subplot(3,1,1);plot(t(2:2:end,1),'color',[1 .0 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_1$ (mm)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
% y
figure(3);
subplot(3,1,2);plot(t(1:2:end,2),'color',[.0 1 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_2$ (mm)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(4);
subplot(3,1,2);plot(t(2:2:end,2),'color',[.0 1 .0],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_2$ (mm)','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
% z
figure(3);
subplot(3,1,3);plot(t(1:2:end,3),'color',[.0 .0 1],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_3$ (mm)','interpreter','latex','fontsize',20);
xlabel('Node indices','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);
figure(4);
subplot(3,1,3);plot(t(2:2:end,3),'color',[.0 .0 1],'linewidth',2);drawnow; axis tight;grid on;
ylabel('$p_3$ (mm)','interpreter','latex','fontsize',20);
xlabel('Node indices','interpreter','latex','fontsize',20);
set(gca,'fontsize',16);

return