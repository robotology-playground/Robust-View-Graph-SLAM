
function Pmax = test_Emat_codes(p1, p2)

kc = zeros(1,5);
K = eye(3);

pts = add_lens_distortion(p1(1:2,:), kc);
pts = pflat(K*pextend(pts));
buckets = test_bucketing(pts(1:2,:), [150 150]);

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);
if size(p1,1)<3; p1 = pextend(p1); end
if size(p2,1)<3; p2 = pextend(p2); end
for ii = 1:10; % how many trials ?
    [p1_sample, p2_sample, pidx] = sample_point_features(p1, p2, buckets, options);
    Evec = calibrated_fivepoint(p1_sample, p2_sample);
    for iiii = 1:size(Evec, 2); % for all possible solutions
        E = reshape(Evec(:,iiii), 3, 3);
        E = E./E(3,3);% normalisation ?
        [R1, R2, t1, t2] = generate_motion_hypothesis(E, C); % generate motion hypothesis
        [P, inlier, err, sol] = resolve_motion_ambiguity(R1, R2, t1, t2, p1, p2, pixtol); % resolve motion ambiguity
        if sum(maxinlier) < sum(inlier);% number of matches critera
            Pmax = P;
            maxinlier = inlier;
            %minerror = norm(err(inlier));
            %pidxmax = pidx;
        end
    end
end

function [p1_sample, p2_sample, idx] = sample_point_features(p1, p2, buckets)
if size(buckets, 1) > 0;
    randind = get_samples(buckets);
else
    randind = randperm(size(p1,2));
end
idx = randind(1:5);
p1_sample = p1(:, idx);
p2_sample = p2(:, idx);

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