
function [P, LS] = test_estimate_extrinsics(p1, p2, M, K1, K2)

% normalise points
p1 = normalise_points(p1, K1);
p2 = normalise_points(p2, K2);

% triangulate using kinematics
%P1 = [eye(3) zeros(3,1)]; % reference camera
%P2 = M(1:3, 1:4);
% xf = intsec2views_midpoint(P1, P2, p1, p2); 
% %posDepth = sum([P{1}(3,:)*U P{2}(3,:)*UU] > 0);

xr = M(1:3, 4); Rr = M(1:3, 1:3);
[r1, r2, D, vis, xf] = test_triangulate(xr, Rr, p1, p2);

% refinement
xc = [xr; rodrigues(Rr)];
xf = xf(1:3, :);
pix = p2(1:2, :);
[xc, xf, vis, x_cov, r, r_cov] = solve_least_squares(pix, xc, xf, eye(3), 200, 1);

% output pairwise geometry structure
P{1} = eye(3,4);
P{2} = [rodrigues(xc(4:6)) xc(1:3)];

% LS output structure
LS.xc = xc;
LS.xf = xf;
LS.x_cov = x_cov;
LS.r = r;
LS.r_cov = r_cov;

% show features
%figure(1); clf;
%line([p1(1, :); p2(1, :)], [p1(2, :); p2(2, :)]); axis equal; drawnow;

% show map
%figure(2); clf;
%plot3(-xf(1,:), xf(3,:), -xf(2,:), '.','color',rand(1,3));%, 'markersize', 2);
%axis equal; grid on; view(3); axis equal; drawnow;

% show depth
%depth = P{2}(3,:)*U;
%figure(3); clf;
%plot(depth);
