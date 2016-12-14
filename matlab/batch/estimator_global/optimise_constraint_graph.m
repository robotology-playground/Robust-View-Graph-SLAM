function [y, Y, sw, x] = optimise_constraint_graph(Ct, C, sw, xs, options)
%
N = length(xs)/6;
%edges=vertcat(C.edge);N=max(edges(:));
if ~isempty(Ct) % Generate trusted information
    Ct = generate_constraint_info_pose(Ct, xs);
end
[y, Y] = initialise_info_matrix(Ct, N);
%
%disp('  ')
disp(['Trusted constraints: ' num2str(length(Ct))])
disp(['ON constraints: ' num2str(sum(sw==1))])
disp(['OFF constraints: ' num2str(sum(sw==0))])
%
C = generate_constraint_info_pose(C, xs); % Initialise estimate with switched ON constraints
[yon, Yon] = generate_joint_info_matrix(C(sw==1), N);
y = y + yon;
Y = Y + Yon;
%
% Apply remaining constraints with residual switching
if options.verbose > 1, figure; end
[y, Y, sw] = constraint_graph_add(y, Y, C, sw, xs, options);
while any(sw)
    [y, Y, sw, ~, converged] = constraint_graph_subtract(y, Y, C, Ct, sw, xs, options);
    if converged, break, end
end
x = recover_moments(y, Y);
%x = Y\y;
%
function [y, Y] = initialise_info_matrix(C, N)
[y, Y] = generate_joint_info_matrix(C, N);
Y(1:6,1:6) = eye(6) * 1e10;