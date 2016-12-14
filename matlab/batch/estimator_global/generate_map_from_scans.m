function generate_map_from_scans(C, pose, scan, kpts, options)

GENERATE_NEW_SCAN = 0;
BUNDLE_ADJUSTMENT = 1;
MAKE_VIDEO = 0;
SCALE = 1000; % (1000-meters) and (1-millimeters)

sigma_r = options.sigma_r/options.K1(1,1); % pixels/focal_length
Rn = (sigma_r.^2)*speye(4); % propagate uncertainty of 1 pixel ?

N = size(pose, 2);
pwg = C_to_pwg(C, N);
for i = 1:N/2
    
    e1 = 2*i-1;
    e2 = 2*i;
    
    % relative motion
    relpose = transform_to_relative_w(pose(:, e2), pose(:, e1));
    R = w2R(relpose(4:6));
    t = relpose(1:3);
    
    % scans in global reference
    if GENERATE_NEW_SCAN == 0
        
        % get local scan and uncertainty
        M = length(scan{e1}.xf)/3;
        xf = reshape(scan{e1}.xf, 3, M);
        P = recover_second_moment(scan{e1}.Y,1);
        % get local scans with acceptable uncertainty
        good = zeros(1, M);
        for j = 1:M
            idx = getindex3(j);
            cov = full(P(idx,idx));
            if trace(cov) < 20^2 %mm^2
                good(j) = 1;
            end
        end
        % move to global reference
        xf = transform_to_global_w(xf(:, good == 1), pose(:, e1));
        
    else
        
        if isempty(pwg{e1, e2}); continue; end;
        p1 = kpts{e1}(pwg{e1, e2}.ind1, :)';
        p2 = kpts{e2}(pwg{e1, e2}.ind2, :)';
        % normalise points
        [K1, K2, kc1, kc2] = get_intrinsics(e1, e2, options);
        % remove lense distortion
        p1 = remove_lens_distortion(p1', kc1, K1); p1 = p1';
        p2 = remove_lens_distortion(p2', kc2, K2); p2 = p2';
        % calibrate
        p1 = K1 \ pextend(p1(1:2, :));
        p2 = K2 \ pextend(p2(1:2, :));
        % generate new local scan
        xf = test_triangulate(R, t, p1, p2);
        % get local scans with acceptable uncertainty
        M = size(xf, 2);
        good = zeros(1, M);
        for j = 1:M
            Hs = test_triangulate_jacobian(R, t, p1(:,j), p2(:,j));
            cov = full(Hs*Rn*Hs');
            if trace(cov) < 500^2 %mm^2
                good(j) = 1;
            end
        end
        % move to global reference
        xf = transform_to_global_w(xf(:, good == 1), pose(:, e1));
        % bundle-adjustment
        if BUNDLE_ADJUSTMENT == 1
            xf = pextend(xf);
            u.pointnr = sum(good);
            u.points{1} = p1(:, good == 1); u.index{1} = 1:u.pointnr;
            u.points{2} = p2(:, good == 1); u.index{2} = 1:u.pointnr;
            xf = modbundle_sparse(xf, {eye(3, 4), [R, t]}, u, 20, 0.001);
        end
        
    end
    
    % plot
    xf = xf / SCALE;
    plot(xf(1,:), xf(3,:), '+');
    xlabel('X (meters)');
    ylabel('Z (meters)');
    grid on, box on, axis equal,
    axis([-1500,1500,0,1500]/SCALE);
    hold on; drawnow,
    %pause(.2)
    
    % video
    if MAKE_VIDEO == 1
        M(i) = getframe(gcf);
    end
    
end

if MAKE_VIDEO == 1
    rate = options.freq/options.steps;
    movie2avi(M, 'aligned_scans.avi', 'fps', rate, 'compression', 'None');
end


%
%
function idx = getindex3(i)
idx = [3*(i-1)+1; 3*(i-1)+2; 3*(i-1)+3];
idx = idx(:);