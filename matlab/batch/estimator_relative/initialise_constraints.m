function [C, Ct, sw] = initialise_constraints(p, options)

% configurations
switch_config;

% generate constraints
sigma_r = options.sigma_r/options.K1(1,1); % pixels/focal_length
R = eye(2)*sigma_r*sigma_r;
ncams = length(p); % ncams = size(p, 1);
npts = length(p(1).s); % npts = size(p{1}, 2);
C = struct('cam', [], 'kpt', [], 'p1', [], 'z', [], 'R', []); % Empty structure fields
k = 0;
for i = 1 : ncams
    
    % supress measurements from images with few matches ?
    if (sum(p(i).s) < options.mincorrnr) %if sum(p{i}(3,:)) < options.mincorrnr
        disp('Image measurements were excluded from constraints at Initialisation');
        disp('This should NOT happen, since poor tracks were removed during Image Processing');
        continue
    end
    
    % initialise measurements with enuogh matches
    for j = 1 : npts
        %if p{i}(3,j) == 1 % status == 1 ?
        if p(i).s(j) == 1 % status == 1 ?
            k = k+1;
            %C(k).cam = i; % cam id (images are numbered 1 : ncams)
            %C(k).kpt = j; % kpt id (points are numbered 1 : npts)
            %C(k).p1 = [p(1).x(j); p(1).y(j)];% p{1}(1:2,j) ;
            %C(k).z = [p(i).x(j); p(i).y(j)];% p{i}(1:2,j);
            %C(k).R = R;
            % the push_back command of matlab
            C(k) = struct( 'cam', i, 'kpt', j, 'p1', [p(1).x(j); p(1).y(j)], ... 
                'z', [p(i).x(j); p(i).y(j)], 'R', R );
        end
    end
    
end

sw = zeros(1, length(C));
%sw = ones(1, length(C));

if TRUST_INIT_LINKS
    Ct = C(sw == 1);
    C = C(sw == 0);
    sw = sw(sw == 0);
else
    Ct = [];
end