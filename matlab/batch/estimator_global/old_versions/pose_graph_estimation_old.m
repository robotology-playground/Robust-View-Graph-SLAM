
tic;
edges=vertcat(C.edge)';
N=max(max(edges));%Number of cameras or images or nodes in view graph

% Formation of A matrix.
m=size(edges,2);
i=[(1:m);(1:m)];i=i(:);
j=edges(:);
s=repmat([-1;1],[m,1]);
k=(j~=1);
%A=sparse(i(k),j(k)-1,s(k),m,N-1);
i=repmat(i(k),1,6);i=i(:); 
j=repmat(j(k),1,6);j=j(:);
s=repmat(s(k),1,6);s=s(:);
block=1; cols = [     6*(j-1)-6+block]; rows = [     6*i-6+block]; data = [     s];
block=2; cols = [cols;6*(j-1)-6+block]; rows = [rows;6*i-6+block]; data = [data;s];
block=3; cols = [cols;6*(j-1)-6+block]; rows = [rows;6*i-6+block]; data = [data;s];
block=4; cols = [cols;6*(j-1)-6+block]; rows = [rows;6*i-6+block]; data = [data;s];
block=5; cols = [cols;6*(j-1)-6+block]; rows = [rows;6*i-6+block]; data = [data;s];
block=6; cols = [cols;6*(j-1)-6+block]; rows = [rows;6*i-6+block]; data = [data;s];
A = sparse(rows, cols, data, 6*m, 6*(N-1));

% measurements information matrix
yy   =1000*repmat([ones(1,3),ones(1,3)],m,1); data = yy';
YY   =sparse(1:6*m,1:6*m,data(:),6*m,6*m);
data =repmat([ones(1,3),ones(1,3)],N-1,1);
Y    =1000*sparse(1:6*(N-1),1:6*(N-1),data(:),6*(N-1),6*(N-1));

% loop this
figure;
sw=false(1,m);
for itr=1:10;%:100;
    
    off=find(sw==0);
    if sum(sw)==m; break; end;
    
    % Residuals
    w           =compute_motion_residuals(q,t,QQ,tt,I);
    s2          =sqrt(sum(w(:,2:4).*w(:,2:4),2)); % norm
    w(:,1)      =2*atan2(s2,w(:,1));
    w(:,1)      =pi_to_pi(w(:,1));
    w(:,2:4)    =w(:,2:4).*repmat(w(:,1)./s2,[1,3]);
    B           =zeros(size(w,1),6);
    B(:,1:3)    =w(:,2:4);
    B(:,4:6)    =w(:,5:7);
    B(isnan(B)) =0;% This tackles the devide by zero problem.
    b           =B';
    b           =b(:);
    
    % residuals switch
    g=zeros(1,sum(~sw));
    for i=1:sum(~sw);
        g(i)=B(off(i),:)*diag(yy(off(i),:))*B(off(i),:)';end
    on=off(g<12.6);
    
    if ~size(on,2);break;end;
    
    plot(g);drawnow;
    fprintf('number of inliers = %d\n',length(on));
    
    % add information
    y=0;
    k=repmat(6*(on-1),6,1)+repmat((1:6)',1,length(on));
    y=y+A(k(:),:)'*YY(k(:),k(:))*b(k(:));
    Y=Y+A(k(:),:)'*YY(k(:),k(:))*A(k(:),:);
    sw(on)=true;
    
    % update state
    x           =Y\y;
    W           =zeros(N,7);
    W(1,:)      =[1 0 0 0 0 0 0];
    W(2:end,2:7)=reshape(x,N-1,6);
    theta       =sqrt(sum(W(:,2:4).*W(:,2:4),2));
    W(:,1)      =cos(theta/2);
    W(:,2:4)    =W(:,2:4).*repmat(sin(theta/2)./theta,[1,3]);
    W(isnan(W)) =0;
    q           =update_rotation(q,W);
    t           =update_position(q,t,W);
end

toc
