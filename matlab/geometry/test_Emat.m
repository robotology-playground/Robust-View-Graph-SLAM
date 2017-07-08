
% Essential matrix
% Tariq Abuhashim - August 2014, iCub

function pwg=test_Emat(cpts1,cpts2,on,K1,K2,kc1,kc2,options,tk)

% Remove points at infinity
vis=remove_points_at_infinity(cpts1(:,on),cpts2(:,on),options.mindisp);
on=find(on);
on=on(vis);
p1=cpts1(:,on);
p2=cpts2(:,on);

% start
pixtol=options.RANSAC_pixtol/K1(1,1);
maxinlier=false;
minerror=options.RANSAC_pixtol; % minimum error at each iteration
pwg=[];

start=tic;

if sum(vis)>options.mincorrnr;  %  minimum number of inliers to trust two-view geometry
    
    % normalise points
    p1=K1\pextend(p1(1:2,:));
    p2=K2\pextend(p2(1:2,:));
    
    % bucketing
    if options.bucketsize(1);
        pts=add_lens_distortion(p1(1:2,:),kc1);
        pts=pflat(K1*pextend(pts));
        buckets=test_bucketing(pts(1:2,:),options.bucketsize(1),options.bucketsize(2));
    end
    
    for ii=1:options.ransac; % how many trials ?
        
        if options.bucketsize(1);
            k=0;buck=1;tri=0;
            ind=zeros(1,size(buckets,1));
            %for buck=1:size(buckets,1)
            while sum(~ind)&&tri<10*size(buckets,1)
                tri=tri+1;
                if buck>size(buckets,1); % checks if end of buckets is reached.
                    if k>4;break;end % finished all buckets and have at least 5 points? Then, exit.
                    buck=2; % finished with less than 5 points? Then, repeat bucketing.
                end;
                data=buckets{buck};
                buck=buck+2; % jump over a bucket, will be revised next step
                if length(data)<3;continue;end; % bucket is empty? Then, move to next one.
                flag=1;idx=0;
                while flag&&idx<10; % if the point is added before, then try at least 10 more times
                    temp=randperm(length(data));
                    i=data(temp(1));
                    flag=sum(ismember(ind,i)); % check if point index is new, flag should be 0
                    idx=idx+1;
                end;
                k=k+1;
                ind(k)=i;
            end
            randind=ind(randperm(k));
        else
            randind=randperm(sum(vis)); % randomise 5 points and best geometry
        end
        p1_sample=p1(:,randind(1:5));
        p2_sample=p2(:,randind(1:5));
        
        % calibrated case
        Evec=calibrated_fivepoint(p2_sample,p1_sample);
        %[E,matches]=test_Fmat_2(p1,p2,K1,K2,opt);Evec=E(:);
        
        for iiii=1:size(Evec,2); % for all possible solutions
            
            E=reshape(Evec(:,iiii),3,3);
            
            [U,~,V]=svd(E);
            if det(U*V')<0;V=-V;end;
            %W = [0 1 0; -1 0 0 ; 0 0 1]; % Olsson
            %Z = [0 -1 0; 1 0 0; 0 0 0];  % Tariq
            W=[0 -1 0;1 0 0 ;0 0 1];  % Harley's book
            %Z = [0 1 0; -1 0 0; 0 0 0];  % Harley's book
            
            P1=[eye(3),zeros(3,1)]; % reference camera
            
            % translations
            %if nargin>8;
            %    t=tk;
            %else
                t=U(:,3);
                %T = U*Z*U'; t1 = [T(3, 2); T(1, 3); T(2, 1)];
                %t = t/norm(t);
                %t = scale.*t;
            %end
            
            
            % first rotation
            R=U*W*V';
            
            % hypothesis #1 : positive translation
            P2=[R,t];
            UU=intsec2views_midpoint(P1,P2,p1,p2);
            P=P2;
            posDepth=sum([P1(3,:)*UU P2(3,:)*UU]>0);
            %xf=test_triangulate(R',t,p1,p2);
            %posDepth=sum(xf(3,:)>0);

            % hypothesis #2 : negative translation
            P2=[R,-t];
            %xf=test_triangulate(R',-t,p1,p2);
            UU=intsec2views_midpoint(P1,P2,p1,p2);
            if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            %if sum(xf(3,:)>0)>posDepth;
                P=P2;
                posDepth=sum([P1(3,:)*UU P2(3,:)*UU]>0);
                %posDepth=sum(xf(3,:)>0);
            end
            
            % second rotation
            R=U*W'*V';
            
            % hypothesis #3 : positive translation
            P2=[R,t];
            %xf=test_triangulate(R',t,p1,p2);
            UU=intsec2views_midpoint(P1,P2,p1,p2);
            if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            %if sum(xf(3,:)>0)>posDepth;
                P=P2;
                posDepth=sum([P1(3,:)*UU P2(3,:)*UU]>0);
                %posDepth=sum(xf(3,:)>0);
            end
            
            % hypothesis #4 : negative translation
            P2=[R,-t];
            %xf=test_triangulate(R',-t,p1,p2);
            UU=intsec2views_midpoint(P1,P2,p1,p2);
            if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            %if sum(xf(3,:)>0)>posDepth;
                P=P2;
                posDepth=sum([P1(3,:)*UU P2(3,:)*UU]>0);
                %posDepth=sum(xf(3,:)>0);
            end
            
            if posDepth %exist('P2', 'var') && ~isempty(P2)
                
                % triangulate all points
                P2=P;
                U=intsec2views_midpoint(P1,P2,p1,p2);
                %R=P2(1:3,1:3);
                %t=P2(1:3,4);
                %xf=test_triangulate(R',t,p1,p2);
                
                % check number of inliers
                err=sqrt(sum((p1-pflat(P1*U)).^2)+sum((p2-pflat(P2*U)).^2));
                mindepth=min(P1(3,:)*U,P2(3,:)*U);
                %[zhat,H]=observation_model(xc,xf);
                %r=zhat+H*[xf(:);xc(:)];
                %err=sum((p2(1:2,:)-reshape(r,2,size(zhat,1)/2)).^2);
                %mindepth=min(xf(3,:));
                inlier=(err<pixtol)&(mindepth>0) ;
                
                %if norm(err(inlier))<minerror;% error exit critera
                if sum(maxinlier)<sum(inlier);% number of matches exit critera
                    %if (norm(err(inlier))<minerror) && (sum(maxinlier)<sum(inlier)) % both
                    Pmax=P2;
                    maxinlier=inlier;
                    minerror=norm(err(inlier));
                end
            end
            
        end
    end
    
    if sum(maxinlier)>options.mincorrnr;
        
        % optimal two-views triangulation
        P2=Pmax;
        U=intsec2views(P1,P2,p1,p2); % optimal two view triangulation
        %U=intsec2views_midpoint(P1,P2,p1,p2);
        %R=P2(1:3,1:3);
        %t=P2(1:3,4);
        %xf=test_triangulate(R',t,p1,p2);
        %clf;plot(xf(1,:),xf(3,:),'+');drawnow;pause;
        
        % pick points based on error and depth sign
        err=sqrt(sum((p1-pflat(P1*U)).^2)+sum((p2-pflat(P2*U)).^2));
        mindepth=min(P1(3,:)*U,P2(3,:)*U);
        %xc=[t;R2w(R')'];
        %[zhat,H]=observation_model(xc,xf);
        %r=zhat+H*[xf(:);xc(:)];
        %err=sum((p2(1:2,:)-reshape(r,2,size(zhat,1)/2)).^2);
        %mindepth=min(xf(3,:));
        maxinlier=err<pixtol&mindepth>0;
        
        % output pairwise geometry structure
        pwg.P=P2; % comment this out is least squares is used
        pwg.inliers=on;
        pwg.e=minerror;
        %pwg.maxinlier = maxinlier;
        pwg.maxinlier=on(maxinlier);
        pwg.opt.ransac=ii;
        %pwg.U = U(:, maxinlier);
        % the points used to calculate E matrix are p1 = p0(ok(vis));
        % the points that produce the best depth map and motion results are
        % p2 = p1(maxinlier)
        % which is:
        % pairwise_geom.inliers = pairwise_geom.inliers(maxinlier);
        
        % bundle-adjustment?
        if 0
            u2.pointnr=size(p1(:,maxinlier),2);
            u2.points{1}=p1(:,maxinlier);
            u2.index{1}=1:u2.pointnr;
            u2.points{2}=p2(:,maxinlier);
            u2.index{2}=1:u2.pointnr;
            [~,P]=modbundle_sparse(U(:,maxinlier),{eye(3,4),Pmax},u2,20,0.001);
            pwg.P=P{2};
        end
        pwg.time=toc(start);
        % show features
        if 0
            clf;
            line([p1(1,maxinlier);p2(1,maxinlier)],[p1(2,maxinlier); ...
                p2(2,maxinlier)]);axis equal
            title(['all visible input points used in estimating the essential matrix : '...
                num2str(sum(maxinlier))]);
            drawnow;
        end
    end
end
