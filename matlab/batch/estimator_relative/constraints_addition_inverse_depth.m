function [y, Y, sw] = constraints_addition_inverse_depth(y, Y, C, sw, xs, options, ncams)
%[y, Y, sw] = constraints_addition_inverse_depth(y, Y, C, sw, xs, options, ncams)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% compute the innovations of the constraints in C that are OFF
off = find(sw == 0);
if length(off)<1; return; end

% recover moments
[x, P] = recover_moments(y, Y);
%r = compute_gate_inverse_depth(x, P, C(off), xs, ncams); %gate of OFF constraints - MATLAB
r = mex_compute_gate_inverse_depth_Mviews(x, P, C(off), xs, ncams); %gate of OFF constraints - cpp
on = off(r < options.gateinnov & r >= 0 & ~isnan(r)); %turn acceptable constraints ON
sw(on) = 1; % Turn ON gated constraints

% insert constraints that passed innovation test
npts = length(xs) - 6*ncams;
%[yon, Yon] = update_info_matrix_inverse_depth(C(on), npts, ncams);
[yon, Yon] = mex_update_info_matrix_Mviews(C(on), npts, ncams);
y = y + yon;
Y = Y + Yon;

% verbose
if options.verbose > 0;
    disp(['Switched ON: ' num2str(length(on)) ', Total: ' ...
        num2str(sum(sw==1)) ' on, ' num2str(sum(sw==0)) ' off']);
end
if options.verbose > 1;
    h1 = subplot(1,2,1); cla(h1);
    plot(r);
    xlabel('Measurements');
    ylabel('Normalised innovations squared')
    title('NIS');
    drawnow; axis tight;
end
