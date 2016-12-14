function [y, Y, sw, x, converged] = constraints_removal_inverse_depth(y, Y, C, Ct, sw, xs, options, ncams)
%[y, Y, sw, x, converged] = constraints_removal_inverse_depth(y, Y, C, Ct, sw, xs, options, ncams)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% compute the innovations of the constraints in C that are ON
on = find(sw == 1);
if length(on)<1; converged = false; x = xs; return; end

% recover moments
%Con = C(on); idx = 1:ncams*6;
%for j = 1:length(Con)
%    idx = [idx, Con(j).kpt+ncams*6];
%end
%idx = sort(unique(idx));
idx = 1:length(xs);
[x, P] = recover_moments(y, Y, idx, ncams);

% print to file
%for i = 1:6*ncams; fprintf(options.fid, '%.4f ', x(i)); end
%fprintf(options.fid, [num2str(length(sw)) ' ' num2str(sum(sw)) '\n']);

%r = compute_gate_inverse_depth(x, P, C(on), xs);%gate of ON constraints - MATLAB
r = mex_compute_gate_inverse_depth_Mviews(x, P, C(on), xs, ncams);%gate of ON constraints - cpp
rmax = max(r(isfinite(r)));
gate = rmax*options.gateratio;

% determine whether to limit minimum gate threshold (based on Ct)
if gate < options.gateresid
    rmaxt = 0;
    if ~isempty(Ct)
        rt = compute_gate_inverse_depth(x, P, Ct, xs);
        %rt = mex_compute_gate_inverse_depth(x, P, Ct, xs);
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

% flag convergence when all ON constraints are less than gate
converged = false;
if rmax < gate
    converged = true;
    return
end
off = on(r > gate | r < 0 | isnan(r));
sw(off) = 0; % Turn OFF gated constraints

% remove constraints that failed residuals test
npts = length(xs) - 6*ncams;
%[yoff, Yoff] = update_info_matrix_inverse_depth(C(off), npts, ncams);
[yoff, Yoff] = mex_update_info_matrix_Mviews(C(off), npts, ncams);
y = y - yoff;
Y = Y - Yoff;

% verbose
if options.verbose ~= 0
    disp(['Switched OFF: ' num2str(length(off)) ', Totals: ' ...
        num2str(sum(sw == 1)) ' on, ' num2str(sum(sw == 0)) ' off']);
end
if options.verbose > 1
    h2 = subplot(1,2,2); cla(h2);
    plot(r);
    xlabel('Measurements');
    ylabel('Normalised residuals squared');
    title('NRS');
    drawnow; axis tight;
end
