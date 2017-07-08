
% Essential matrix
% Tariq Abuhashim - August 2014, iCub

function pwg = test_Emat_v2(p1, p2, C, options)


pwg = [];
if size(p1,2) < options.mincorrnr;  return; end;%  minimum number of inliers to trust two-view geometry

start = tic;
% start
[K1, K2, kc1, kc2] = get_intrinsics(options, C.edge(1), C.edge(2));
pixtol = options.RANSAC_pixtol/K1(1, 1);
maxinlier = false;
minerror = options.RANSAC_pixtol; % minimum error at each iteration

% bucketing
if options.bucketsize(1);
    pts = add_lens_distortion(p1(1:2,:),kc1);
    pts = pflat(K1*pextend(pts));
    buckets = test_bucketing(pts(1:2,:), options.bucketsize(1), options.bucketsize(2));
end

for ii = 1:options.ransac; % how many trials ?
    
    % randomise 5 points and best geometry
    if options.bucketsize(1);
        randind = get_samples(buckets, options);
    else
        randind = randperm(sum(vis));
    end
    p1_sample = p1(:, randind(1:5));
    p2_sample = p2(:, randind(1:5));
    
    % calibrated case
    Evec = calibrated_fivepoint(p2_sample, p1_sample);
    
    for iiii = 1:size(Evec,2); % for all possible solutions
        
        E = reshape(Evec(:,iiii),3,3);
        [U,~,V] = svd(E);
        if det(U*V')<0; V=-V; end;
        W = [0 -1 0; 1 0 0 ; 0 0 1];  % Harley's book
        
        % reference camera
        P1 = [eye(3), zeros(3, 1)];
        
        % translations
        t = C.t';
        %t = U(:,3);
        
        % first rotation
        R = U*W*V';
        
        % hypothesis #1 : positive translation
        P2 = [R, t];
        UU = intsec2views_midpoint(P1, P2, p1, p2);
        P = P2;
        posDepth = sum([P1(3,:)*UU P2(3,:)*UU]>0);
        
        % hypothesis #2 : negative  translation
        P2 = [R, -t];
        UU = intsec2views_midpoint(P1, P2, p1, p2);
        if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            P = P2;
            posDepth = sum([P1(3,:)*UU P2(3,:)*UU]>0);
        end
        
        % second rotation
        R = U*W'*V';
        
        % hypothesis #3 : positive translation
        P2 = [R, t];
        UU = intsec2views_midpoint(P1, P2, p1, p2);
        if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            P = P2;
            posDepth = sum([P1(3,:)*UU P2(3,:)*UU]>0);
        end
        
        % hypothesis #4 : negative translation
        P2 = [R, -t];
        UU = intsec2views_midpoint(P1, P2, p1, p2);
        if sum([P1(3,:)*UU P2(3,:)*UU]>0)>posDepth
            P = P2;
            posDepth = sum([P1(3,:)*UU P2(3,:)*UU]>0);
        end
        
        if posDepth %exist('P2', 'var') && ~isempty(P2)
            % triangulate all points
            P2 = P;
            U = intsec2views_midpoint(P1, P2, p1, p2);
            % check number of inliers
            err = sqrt(sum((p1-pflat(P1*U)).^2) + sum((p2-pflat(P2*U)).^2));
            mindepth = min(P1(3,:)*U,P2(3,:)*U);
            inlier = (err<pixtol)&(mindepth>0) ;
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
    P2 = Pmax;
    U = intsec2views(P1, P2, p1, p2); % optimal two view triangulation
    
    % pick points based on error and depth sign
    err = sqrt(sum((p1-pflat(P1*U)).^2)+sum((p2-pflat(P2*U)).^2));
    mindepth = min(P1(3,:)*U,P2(3,:)*U);
    maxinlier = err<pixtol&mindepth>0;
    
    % output pairwise geometry structure
    pwg.P = P2; % comment this out is least squares is used
    %pwg.inliers = on;
    pwg.e = minerror;
    %pwg.maxinlier = maxinlier;
    pwg.maxinlier = maxinlier;
    pwg.opt.ransac = ii;
    
    % bundle-adjustment?
    if 0
        u2.pointnr = size(p1(:, maxinlier), 2);
        u2.points{1} = p1(:, maxinlier);
        u2.index{1} = 1:u2.pointnr;
        u2.points{2} = p2(:, maxinlier);
        u2.index{2} = 1:u2.pointnr;
        [~, P] = modbundle_sparse(U(:, maxinlier),{eye(3, 4),Pmax},u2, 20, 0.001);
        pwg.P = P{2};
    end
    
    pwg.time = toc(start);
    
    % show features
    if 0
        clf;
        line([p1(1, maxinlier); p2(1, maxinlier)],...
            [p1(2, maxinlier); p2(2, maxinlier)]);
        axis equal
        title(['all visible input points used in estimating the essential matrix : '...
            num2str(sum(maxinlier))]);
        drawnow;
    end
end

function randind = get_samples(buckets, options)
if options.bucketsize(1);
    k = 0;
    buck = 1;
    tri = 0;
    ind = zeros(1,size(buckets, 1));
    while sum(~ind) && tri < 10*size(buckets, 1)
        tri = tri + 1;
        if buck > size(buckets, 1); % checks if end of buckets is reached.
            if k > 4;
                break; % finished all buckets and have at least 5 points? Then, exit.
            end
            buck = 2; % finished with less than 5 points? Then, repeat bucketing.
        end
        data = buckets{buck};
        buck = buck + 2; % jump over a bucket, will be revised next step
        if length(data) < 3;
            continue;
        end % bucket is empty? Then, move to next one.
        flag = 1;
        idx = 0;
        while flag && idx < 10; % if the point is added before, then try at least 10 more times
            temp = randperm(length(data));
            i = data(temp(1));
            flag = sum(ismember(ind, i)); % check if point index is new, flag should be 0
            idx = idx + 1;
        end
        k = k + 1;
        ind(k) = i;
    end
    randind = ind(randperm(k));
end
return