

% Essential matrix
% Tariq Abuhashim - August 2014, iCub

function pwg = test_Mransac(cpts1, cpts2, ok, K1, K2, mindisp, RANSAC_pixtol, mincorrnr, M)

% Remove points at infinity
vis = remove_points_at_infinity(cpts1(:, ok), cpts2(:, ok), mindisp);
ok = find(ok);
p1 = cpts1(:, ok(vis));
p2 = cpts2(:, ok(vis));

% start
pixtol = RANSAC_pixtol/K1(1,1);
maxinlier = false;
minerror = RANSAC_pixtol; % minimum error at each iteration
pwg = [];
s_variance = [[.01;.01;.01];[.01;.01;.01].*pi/180]; % kinematics uncertainty
if sum(vis) > mincorrnr;  %  minimum number of inliers to trust two-view geometry
    
    % normalise points
    p1 = normalise_points(p1, K1);
    p2 = normalise_points(p2, K2);
    
    % sample motion prior
    n_samples = 2000;
    r = rodrigues(M(1:3, 1:3));
    t = M(1:3, 4);
    x = normrnd(repmat([t;r],1,n_samples),repmat(s_variance,1,n_samples),[6 n_samples]);
    %x = [t;r];
    for ii = 1:size(x,2) % how many trials ?
        
        P{1} = [eye(3) zeros(3,1)]; % reference camera
        t = x(1:3,ii);
        R = rodrigues(x(4:6,ii));
        
        % hypothesis #1 : positive translation
        P{2} = [R t];
        UU = intsec2views_midpoint(P{1}, P{2}, p1, p2);
        posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        P2 = P{2};
        
        % hypothesis #2 : negative translation
        P{2} = [R -t];
        UU = intsec2views_midpoint(P{1}, P{2}, p1, p2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2};
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        % hypothesis #3 : positive translation
        P{2} = [R' t];
        UU = intsec2views_midpoint(P{1}, P{2}, p1, p2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2};
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        % hypothesis #4 : negative translation
        P{2} = [R' -t];
        UU = intsec2views_midpoint(P{1}, P{2}, p1, p2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2};
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        if posDepth %exist('P2', 'var') && ~isempty(P2)
            
            % triangulate all points
            P{1} = [eye(3) zeros(3,1)];
            P{2} = P2;
            U = intsec2views_midpoint(P{1}, P{2}, p1, p2);
            
            % check number of inliers
            err = sqrt( sum((p1-pflat(P{1}*U)).^2) + sum((p2-pflat(P{2}*U)).^2) );
            mindepth = min(P{1}(3,:)*U,P{2}(3,:)*U);
            inlier = ( err < pixtol ) & ( mindepth > 0 ) ;
            
            %if norm(err(inlier))<minerror;% error exit critera
            if sum(maxinlier)<sum(inlier);% number of matches exit critera
            %if (norm(err(inlier))<minerror) && (sum(maxinlier)<sum(inlier)) % both
                Pmax = P;
                maxinlier = inlier;
                minerror = norm(err(inlier));
            end
        end
    end
end

if sum(maxinlier) > mincorrnr;
    
    % optimal two-views triangulation
    P = Pmax;
    %U = intsec2views(P{1},P{2},p1,p2); % optimal two view triangulation
    U = intsec2views_midpoint(P{1},P{2},p1,p2);
    
    % pick points based on error and depth sign
    err = sqrt( sum((p1-pflat(P{1}*U)).^2) + sum((p2-pflat(P{2}*U)).^2) );
    mindepth = min(P{1}(3,:)*U,P{2}(3,:)*U);
    maxinlier = err < pixtol & mindepth > 0;
    u.pointnr = size(p1(:, maxinlier), 2);
    u.points{1} = p1(:, maxinlier);
    u.index{1} = 1:u.pointnr;
    u.points{2} = p2(:, maxinlier);
    u.index{2} = 1:u.pointnr;
    U = U(:, maxinlier);
    
    % refinement
    %[U, P] = modbundle_sparse(U, P, u, 20, .1);
    x0 = [P{2}(:,4); rodrigues(P{2}(1:3,1:3))];
    pix = p2(1:2,maxinlier);
    % first with gate = 1
    [xc, U, vis2] = solve_least_squares(pix, x0, U(1:3,:), eye(3), 1, 1);
    % then Uith gate = 0
    [xc, U] = solve_least_squares(pix(:, vis2), xc, U, eye(3), 200, 0);
    %[xc, U] = solve_least_squares(pix, x0, U(1:3,:), eye(3), 200, 0);

    P{2} = [rodrigues(xc(4:6)) xc(1:3)];
    
    % output pairwise geometry structure
    pwg.P = P;
    pwg.U = U;
    pwg.inliers = ok(vis);
    pwg.e = minerror;
    pwg.maxinlier = maxinlier;
    
    % show features
    %figure(1); clf; 
    %line([p1(1, maxinlier); p2(1, maxinlier)], [p1(2, maxinlier); ...
    %    p2(2, maxinlier)]); axis equal
    %title(['all visible input points used in estimating the essential matrix : '...
    %    num2str(sum(maxinlier))]); drawnow;
    
    % show map
    U = intsec2views_midpoint(P{1},P{2},p1(:,maxinlier), p2(:,maxinlier)); 
    figure(2); clf;  
    plot3(-U(1,:), U(3,:), -U(2,:), '.','color',rand(1,3));%, 'markersize', 2);
    axis equal; grid on; view(3); axis equal; drawnow;
    
    depth = P{2}(3,:)*U;
    figure(3); clf;
    plot(depth);
    
    pause;
end