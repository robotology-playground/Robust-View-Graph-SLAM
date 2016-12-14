function QQ=quaternion_graph_distance(xs,edges)

N=length(xs)/6;
RR=zeros(3,3,N);
for i = 1:N
    RR(:,:,i)=w2R(xs((i-1)*6+(4:6)));
end
Q=[RR(1,1,:)+RR(2,2,:)+RR(3,3,:)-1, RR(3,2,:)-RR(2,3,:),RR(1,3,:)-RR(3,1,:),RR(2,1,:)-RR(1,2,:)]/2;
Q=reshape(Q,4,size(Q,3),1)';
Q(:,1)=sqrt((Q(:,1)+1)/2);
Q(:,2:4)=(Q(:,2:4)./repmat(Q(:,1),[1,3]))/2;

i=edges(1,:);
j=edges(2,:);
    
% w=Qij*Qi
%w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
%    repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
%    [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

% w=inv(Qj)*w=inv(Qj)*Qij*Qi
QQ=[ (-Q(j,1).*Q(i,1)-sum(Q(j,2:4).*Q(i,2:4),2)),...  %scalar terms
    repmat(-Q(j,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*Q(j,2:4) + ...   %vector terms
    [Q(j,3).*Q(i,4)-Q(j,4).*Q(i,3),Q(j,4).*Q(i,2)-Q(j,2).*Q(i,4),Q(j,2).*Q(i,3)-Q(j,3).*Q(i,2)] ];   %cross product terms

% for i=1:size(Q,1)
%     R=real(q2R(Q(i,:)));
%     R2w(R);
% end