function [C, use, J] = test_Emat_v3(C, p1, p2, use, options)
%[C, J] = test_Emat_v3(C, p1, p2, use, options)
%
% Robust essential matrix estimation from noisy image points using RANSAC
%
% INPUTS :
%   C : constraints
%   p1,p2 : feature matches
%   use : index of which features to use
%   options : options
%
% OUPUTS :
%   C : updated constraints using epipolar geometry
%   J : jacobian of the essential matrix
%
% Few other attempts may be found at:
% 1- The current code benifits from Stew´enius's implementation of Nestor's 5-points algorithm
%       Available at : http://vis.uky.edu/~stewe/FIVEPOINT/
% 2- Nghia Ho has implemented a C++ code, but could potentially contain errors (check the comments)
%       Available at : http://nghiaho.com/?p=1675
% 3- Manolis Lourakis also has a robust version of the fundamental matrix
%       Available at : http://users.ics.forth.gr/~lourakis/fundest/
% 4- Avi Sigh has also great complete tutorial with codes for visual odometry using monocular and
%    stereo
%       Available at : https://avisingh599.github.io/vision/visual-odometry-full/
%                      http://avisingh599.github.io/research/
% 5- Visual inertial odometry at MATLAB file exchange
%       Available at : http://www.mathworks.com/matlabcentral/fileexchange/43218-visual-inertial-odometry
%
% This code utilises features selection through spatial division of the image into multiple regions.
% Points are then collected in such away to improve the sampled geometry.
% This is a simple rather a smart way of doing this.
%
% Tariq Abuhashim, 2015.
%
% iCub - Koroibot
%

s = RandStream('mlfg6331_64','Seed',0);
RandStream.setGlobalStream(s);

p1 = p1(:, use);
p2 = p2(:, use);
[K, kc] = get_intrinsics(options, C.edge(1));
pixtol = options.RANSAC_pixtol/K(1, 1);
maxinlier = false;
%minerror = options.RANSAC_pixtol; % minimum error at each iteration

% fundamental matrix
%[F, idxOutliers] = fundest(p1(1:2, :)', p2(1:2, :)', 0.7, 1, 'sampson_err');
%inliers = ones(1, size(p1, 2));
%inliers(idxOutliers) = 0;
%p1 = p1(:, inliers == 1);
%p2 = p2(:, inliers == 1);
%use = use(inliers == 1);

% bucketing (features selection)
buckets = [];
if options.bucketsize(1) > 0;
    pts = add_lens_distortion(p1(1:2,:), kc);
    pts(1,:)=K(1,1)*pts(1,:)+K(1,3);
    pts(2,:)=K(2,2)*pts(2,:)+K(2,3);
    %pts = pflat(K*pextend(pts));
    %im1 = get_image(C.edge(1), C.edge(2), options);
    buckets = test_bucketing(pts(1:2,:), options.bucketsize);%, im1);
end

if size(p1,1)<3; p1 = pextend(p1); end
if size(p2,1)<3; p2 = pextend(p2); end
for ii = 1:options.ransac; % how many trials ?
    % randomise 5 points and best geometry
    [p1_sample, p2_sample, pidx] = sample_point_features(p1, p2, buckets, options);
                %clf; plot(p1_sample(1,:), p1_sample(2,:), '+');
                %hold on; plot(p2_sample(1,:), p2_sample(2,:), 'ro'); pause;
    % calibrated case
    Evec = calibrated_fivepoint(p1_sample, p2_sample);
    for iiii = 1:size(Evec, 2); % for all possible solutions
        E = reshape(Evec(:,iiii), 3, 3);
        E = E./E(3,3);% normalisation ?
        [R1, R2, t1, t2] = generate_motion_hypothesis(E, C); % generate motion hypothesis
        [P, inlier, err, sol] = resolve_motion_ambiguity(R1, R2, t1, t2, p1, p2, pixtol); % resolve motion ambiguity
        % test the hypothesis
        %if norm(err(inlier)) < minerror;% minimum error critera
        if sum(maxinlier) < sum(inlier);% number of matches critera
        %if (norm(err(inlier)) < minerror) && (sum(maxinlier)<sum(inlier)) % both
            Pmax = P;
            maxinlier = inlier;
            minerror = norm(err(inlier));
            pidxmax = pidx;
            %idx = iiii;
            %solmin = sol;
            %p1min = p1_sample(1:2,:);
            %p2min = p2_sample(1:2,:);
        end
    end
end

% bundle-adjustment ?
if 0
    %xc = [Pmax(1:3,4); R2w(Pmax(1:3,1:3))'];
    %r1 = test_triangulate_inverse_depth(p1, p2, xc);
    %z1 = get_scan_from_range(p1,1./r1); U = pextend(z1);
    P1 = eye(3,4); P2 = [Pmax(1:3,1:3)' -Pmax(1:3,1:3)'*Pmax(1:3,4)];
    U = intsec2views(P1, P2, p1, p2); % optimal two view triangulation
    % pick points based on error and depth sign
    err = sqrt(sum((p1-pflat(P1*U)).^2)+sum((p2-pflat(P2*U)).^2));
    mindepth = min(P1(3,:)*U,P2(3,:)*U);
    maxinlier2 = err<pixtol&mindepth>0;
    u2.pointnr = size(p1(:,maxinlier2), 2);
    u2.points{1} = p1(:,maxinlier2);
    u2.index{1} = 1:u2.pointnr;
    u2.points{2} = p2(:,maxinlier2);
    u2.index{2} = 1:u2.pointnr;
    [~, P] = modbundle_sparse(U(:,maxinlier2),{P1,P2},u2, 20, 0.001);
    Pmax = [P{2}(1:3,1:3)' -P{2}(1:3,1:3)'*P{2}(1:3,4)];
end
        
% output pairwise geometry structure
C.r = Pmax(1:3,1:3);
C.q = w2q(R2w(C.r));
C.t = Pmax(1:3,4)';
C.maxinlier = maxinlier;
C.e = minerror;
C.p_idx = pidxmax;

%use = use(maxinlier == 1);

% get the jacobian of the essential matrix
if nargout>2
    %J = Emat_jacobian(p1min, p2min, idx, solmin);
    %J*(((options.sigma_r/410)^2)*eye(20))*J'
    J = [];
end

function [p1_sample, p2_sample, idx] = sample_point_features(p1, p2, buckets, options)
if size(buckets, 1) > 0;
    randind = get_samples(buckets);
else
    randind = randperm(size(p1,2));
end
idx = randind(1:5);
p1_sample = p1(:, idx);
p2_sample = p2(:, idx);
% plot sampled points
if 0;%options.verbose > 1
    pts = add_lens_distortion(p1_sample, kc1);
    pts = pflat(K1*pextend(pts));
    plot(pts(1,:),pts(2,:),'rO','markersize',10,'linewidth',2);
    drawnow;
end

function randind = get_samples(buckets) % an older version is located at the end of this file
k = 0;
buck = 1; % start with odd buckets first
tri = 0;
ind = zeros(1, max(5, size(buckets, 1))); % because minimum 5 points are needed
while sum(~ind) && tri < min(2000, 10*size(buckets, 1))
    tri = tri + 1;
    if buck > size(buckets, 1); % checks if end of buckets is reached.
        if k > 4;
            break; % finished all buckets and have at least 5 points? Then, exit.
        end
        buck = 1+mod(buck,2); % finished with less than 5 points? Then, repeat bucketing using even buckets
    end
    data = buckets{buck};
    buck = buck + 2; % jump over a bucket, will be revised next step
    if length(data) < 3;
        continue;
    end % bucket is empty? Then, move to next one.
    flag = 1;
    idx = 0;
    temp = randperm(length(data));
    while flag && idx<length(data); % if the point is added before, then try at least 10 more times
        idx = idx + 1;
        i = data(temp(idx));
        %i = data(idx);
        flag = sum(ismember(ind, i)); % check if point index is new, flag should be 0
    end
    k = k + 1;
    ind(k) = i;
end
randind = ind(randperm(k));

function [R1, R2, t1, t2] = generate_motion_hypothesis(E, C)
% re-enforce rank 2 constraint on essential matrix
[U,D,V] = svd(E,0);
E = U*diag([D(1,1) D(2,2) 0])*V';
% SVD decomposition of the essential matrix
[U,~,V] = svd(E,0);
if det(U*V') < 0; V = -V; end;
W = [0 -1 0; +1 0 0 ; 0 0 +1]; % Harley's book
%Z = [0 +1 0; -1 0 0 ; 0 0  0]; % Harley's book
R1 = U*W*V';
R2 = U*W'*V';
t1 = U(:,3);
if nargin>1
    if isfield(C, 't');
        t1 = C.t';
    end
end
t2 = -t1;
% assure determinant to be positive
if det(R1) < 0; R1 = -R1; end;
if det(R2) < 0; R2 = -R2; end;

function [P, inlier, err, sol] = resolve_motion_ambiguity(R1, R2, t1, t2, p1, p2, pixtol)
% case 1: {R1,t1}
xc = [t1; R2w(R1)'];
r1 = test_triangulate_inverse_depth(p1, p2, xc);
%xc = transform_to_relative_w(zeros(6,1), xc); 
%xc = [-R1'*t1; R2w(R1')'];
%r2 = test_triangulate_inverse_depth(p2, p1, xc);
posDepth = sum(r1 > 0);
P = [R1, t1]; sol = 1;
% case 2: {R1,t2}
xc = [t2; R2w(R1)'];
r1 = test_triangulate_inverse_depth(p1, p2, xc);
%xc = transform_to_relative_w(zeros(6,1), xc); 
%xc = [-R1'*t2; R2w(R1')'];
%r2 = test_triangulate_inverse_depth(p2, p1, xc);
if sum(r1 > 0) > posDepth
    P = [R1, t2]; sol = 2;
    posDepth = sum(r1 > 0);
end
% case 3: {R2,t1}
xc = [t1; R2w(R2)'];
r1 = test_triangulate_inverse_depth(p1, p2, xc);
%xc = transform_to_relative_w(zeros(6,1), xc); 
%xc = [-R2'*t1; R2w(R2')'];
%r2 = test_triangulate_inverse_depth(p2, p1, xc);
if sum(r1 > 0) > posDepth
    P = [R2, t1]; sol = 2;
    posDepth = sum(r1 > 0);
end
% case 4: {R2,t2}
xc = [t2; R2w(R2)'];
r1 = test_triangulate_inverse_depth(p1, p2, xc);
%xc = transform_to_relative_w(zeros(6,1), xc); 
%xc = [-R2'*t2; R2w(R2')'];
%r2 = test_triangulate_inverse_depth(p2, p1, xc);
if sum(r1 > 0) > posDepth
    P = [R2, t2]; sol = 2;
    posDepth = sum(r1 > 0);
end
if posDepth
    R = P(1:3,1:3); t = P(1:3,4);
    % minimum depth
    xc = [t; R2w(R)'];
    r1 = test_triangulate_inverse_depth(p1, p2, xc);
    %xc = transform_to_relative_w(zeros(6,1), xc);
    %xc = [-R'*t; R2w(R')'];
    %r2 = test_triangulate_inverse_depth(p2, p1, xc);
    % error
    z1 = get_scan_from_range(p1,1./r1);
    z2 = z1;
    z2(1,:) = z1(1,:) - t(1);
    z2(2,:) = z1(2,:) - t(2);
    z2(3,:) = z1(3,:) - t(3);
    z2 = pflat(R'*z2);
    err = sqrt(sum((p2(1:2,:)-z2(1:2,:)).^2));
    %err = sqrt(sum((p1(1:2,:)-z1(1:2,:)).^2) + sum((p2(1:2,:)-z2(1:2,:)).^2));
    % inliers
    inlier = (err < pixtol) & (r1 > 0);%&(r2 > 0);
end

function p = get_scan_from_range(p,r)
if size(p,1) < 3; p = pextend(p); end;
rim = sqrt(sum(p(1:2,:).^2) + 1);
d = r./rim;
p(1,:) = p(1,:).*d;
p(2,:) = p(2,:).*d;
p(3,:) = p(3,:).*d;

function J = Emat_jacobian(p1, p2, idx, solmin)
J = [];
for i = 1:10;
    Ji = numerical_jacobian_i(@test_Evec,[],i,[],...
        p1(:,1), p1(:,2), p1(:,3), p1(:,4), p1(:,5),...
        p2(:,1), p2(:,2), p2(:,3), p2(:,4), p2(:,5), idx, solmin);
    J = [J Ji];
end

function out = test_Evec(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,idx,solmin)
p1_sample = [p0 p1 p2 p3 p4]; p1_sample = pextend(p1_sample);
p2_sample = [p5 p6 p7 p8 p9]; p2_sample = pextend(p2_sample);
Evec = calibrated_fivepoint(p1_sample, p2_sample);
E = reshape(Evec(:,idx),3,3);
% re-enforce rank 2 constraint on essential matrix
[U,D,V] = svd(E,0);
E = U*diag([D(1,1) D(2,2) 0])*V';
% SVD decomposition of the essential matrix
[U,~,V] = svd(E,0);
if det(U*V') < 0; V = -V; end;
W = [0 -1 0; +1 0 0 ; 0 0 +1]; % Harley's book
%Z = [0 +1 0; -1 0 0 ; 0 0  0]; % Harley's book
R1 = U*W*V';
R2 = U*W'*V';
t = U(:,3);
% assure determinant to be positive
if det(R1) < 0; R1 = -R1; end;
if det(R2) < 0; R2 = -R2; end;
switch solmin
    case 1;
        w = R2w(R1);
    case 2;
        w = R2w(R1); t = -t;
    case 3;
        w = R2w(R2);
    case 4;
        w = R2w(R2); t = -t;
end
%out = Evec(:,idx); % Jacobian of essential matrix
out = [t(:); w(:)]; % Jacobian of motion parameters
        
% function [P, inlier, err, sol] = resolve_motion_ambiguity(R1, R2, t1, t2, p1, p2, pixtol)
% % case 1: {R1,t1}
% xf1 = test_triangulate(R1, t1, p1, p2);
% xf2 = transform_to_relative_w(xf1, [t1; R2w(R1)']);
% posDepth = sum([xf1(3,:) xf2(3,:)] > 0);
% P = [R1, t1]; sol = 1;
% % case 2: {R1,t2}
% xf1 = test_triangulate(R1, t2, p1, p2);
% xf2 = transform_to_relative_w(xf1, [t2; R2w(R1)']);
% if sum([xf1(3,:) xf2(3,:)] > 0) > posDepth
%     P = [R1, t2]; sol = 2;
%     posDepth = sum([xf1(3,:) xf2(3,:)] > 0);
% end
% % case 3: {R2,t1}
% xf1 = test_triangulate(R2, t1, p1, p2);
% xf2 = transform_to_relative_w(xf1, [t1; R2w(R2)']);
% if sum([xf1(3,:) xf2(3,:)] > 0) > posDepth
%     P = [R2, t1]; sol = 3;
%     posDepth = sum([xf1(3,:) xf2(3,:)] > 0);
% end
% % case 4: {R2,t2}
% xf1 = test_triangulate(R2, t2, p1, p2);
% xf2 = transform_to_relative_w(xf1, [t2; R2w(R2)']);
% if sum([xf1(3,:) xf2(3,:)] > 0) > posDepth
%     P = [R2, t2]; sol = 4;
%     posDepth = sum([xf1(3,:) xf2(3,:)] > 0);
% end
% if posDepth
%     % triangulate all points
%     R = P(1:3, 1:3); t = P(1:3, 4);
%     xf1 = test_triangulate(R, t, p1, p2); z1 = pflat(xf1);
%     xf2 = transform_to_relative_w(xf1, [t; R2w(R)']); z2 = pflat(xf2);
%     % check number of inliers
%     %err = sqrt(sum((p1(1:2,:)-z1(1:2,:)).^2) + sum((p2(1:2,:)-z2(1:2,:)).^2));
%     err = sqrt(sum((p2(1:2,:)-z2(1:2,:)).^2));
%     mindepth = sum([xf1(3,:) xf2(3,:)] > 0);
%     inlier = (err < pixtol) & (mindepth > 0);
% end

% function randind = get_samples(buckets, options) % old version of features bucketing
%
% ODD_FIRST = 1; % (1)-odd first, (2)-even first, (0)-first-bucket first
%
% if options.bucketsize(1);
%     k = 0;
%     if ODD_FIRST == 2
%         buck = 2; % start with even buckets first
%     else
%         buck = 1; % start with odd buckets first
%     end
%     tri = 0;
%     ind = zeros(1,size(buckets, 1));
%     while sum(~ind) && tri < 10*size(buckets, 1)
%         tri = tri + 1;
%         if buck > size(buckets, 1); % checks if end of buckets is reached.
%             if k > 4;
%                 break; % finished all buckets and have at least 5 points? Then, exit.
%             end
%             if ODD_FIRST == 1
%                 buck = 2; % finished with less than 5 points? Then, repeat bucketing using even buckets
%             else
%                 buck = 1; % finished with less than 5 points? Then, repeat bucketing from start
%             end
%         end
%         data = buckets{buck};
%
%         if ODD_FIRST > 0
%             buck = buck + 2; % jump over a bucket, will be revised next step
%         else
%             buck = buck + 1; % do not jump over a bucket
%         end
%
%         if length(data) < 3;
%             continue;
%         end % bucket is empty? Then, move to next one.
%         flag = 1;
%         idx = 0;
%         while flag && idx < 10; % if the point is added before, then try at least 10 more times
%             temp = randperm(length(data));
%             i = data(temp(1));
%             flag = sum(ismember(ind, i)); % check if point index is new, flag should be 0
%             idx = idx + 1;
%         end
%         k = k + 1;
%         ind(k) = i;
%     end
%     randind = ind(randperm(k));
% end