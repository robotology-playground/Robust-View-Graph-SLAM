function [y, Y, sw] = constraint_graph_add(y, Y, C, sw, xs, options)
% Insert constraints that pass innovation test
[x, P] = recover_moments(y, Y);
if options.verbose > 1, plot_pose_graph_w(x); end

off = find(sw == 0); % compute innovations (of OFF constraints)
%r = compute_residuals(x, P, C(off), xs);
r = mex_compute_gate_graph(x, P, C(off), xs);
on = off(r < options.gateinnov); % turn acceptable constraints ON
sw(on) = 1;

N = length(xs)/6;
[yon, Yon] = generate_joint_info_matrix(C(on), N);
Y = Y + Yon;    
y = y + yon;

if options.verbose ~= 0
    disp(['Switched ON: ' num2str(length(on)) ', Totals: ' ...
        num2str(sum(sw==1)) ' on, ' num2str(sum(sw==0)) ' off'])
end