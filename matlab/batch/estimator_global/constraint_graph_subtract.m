function [y, Y, sw, x, converged] = constraint_graph_subtract(y, Y, C, Ct, sw, xs, options)
% Turn off constraints based on residuals.
[x,P] = recover_moments(y, Y);
if options.verbose > 1;
    clf; 
    plot_pose_graph_w(x); 
end

% Compute residuals of constraints currently ON
on = find(sw == 1);
%r = compute_residuals(x, P, C(on), xs);
r = mex_compute_gate_graph(x, P, C(on), xs);
rmax = max(r);
gate = rmax * options.gateratio;

% Determine whether to limit minimum gate threshold (based on Ct)
if gate < options.gateresid
    rmaxt = 0;
    if ~isempty(Ct)
        rt = compute_residuals(x, P, Ct, xs);
        rmaxt = max(rt);
        disp(['Max trusted residual: ' num2str(rmaxt)])
    end
    if rmaxt < options.gatetrust
        gate = options.gateresid; % limit gate to minimum threshold
    end
end
if options.verbose ~= 0
    disp(['Max residual: ' num2str(rmax) ', Gate: ' num2str(gate)])
end

% Flag convergence when all ON constraints are less than gate
converged = false;
if rmax < gate
    converged = true;
    return
end

% Turn OFF gated constraints
off = on(r > gate);
sw(off) = 0;
N = length(xs)/6;
[yoff, Yoff] = generate_joint_info_matrix(C(off), N);
Y = Y - Yoff;    
y = y - yoff;
if options.checkrank % FIXME: monitor spanning tree to avoid rank deficiency
end

if options.verbose ~= 0
    disp(['Switched OFF: ' num2str(length(off)) ', Totals: ' ...
        num2str(sum(sw==1)) ' on, ' num2str(sum(sw==0)) ' off'])
end
