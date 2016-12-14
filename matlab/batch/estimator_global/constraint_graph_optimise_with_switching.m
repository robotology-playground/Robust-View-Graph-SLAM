function [y, Y, sw, x] = constraint_graph_optimise_with_switching(y, Y, C, sw, xs, options)
% This function presumes that {y,Y} already contains the information from
% C(sw==1).
if nargin < 7
    options = constraint_graph_optimise_set_options; % default configuration
end
N = length(xs)/6;

for i = 1:options.iterations
    % Insert constraints that pass innovation test
    clear P
    [x,P] = recover_moments(y, Y);
    if options.verbose > 1, plot_pose_graph(x); end
    off = find(sw == 0);
    r = compute_residuals(x, P, C(off), xs);
    on = off(r < options.gateinnov);
    sw(on) = 1;
    [yon, Yon] = generate_joint_info_matrix(C(on), N);
    Y = Y + Yon;    
    y = y + yon;
    if options.verbose ~= 0, disp(['Switched ON: ' num2str(length(on))]), end
    
    % Remove constraints that fail residual test
    while 1
        clear P
        [x,P] = recover_moments(y, Y);
        if options.verbose > 1, plot_pose_graph(x); end
        on = find(sw == 1);
        r = compute_residuals(x, P, C(on), xs);
        rmax = max(r);
        gate = rmax * options.gateratio;
        if gate < options.gateresid, gate = options.gateresid; end
        if options.verbose ~= 0
            disp(['Max residual: ' num2str(rmax) ', Gate: ' num2str(gate)])
        end
        if rmax < options.gateresid, break, end            
        off = on(r > gate);
        sw(off) = 0;
        [yoff, Yoff] = generate_joint_info_matrix(C(off), N);
        Y = Y - Yoff;    
        if options.checkrank % FIXME: monitor spanning tree to avoid rank deficiency
        end
        y = y - yoff;
        if options.verbose ~= 0, disp(['Switched OFF: ' num2str(length(off))]), end
    end
end
