function [x, P, C, sw, GATE_PROB_RESID] = PwgOptimiser(C, Ct, xs, sw, ncams, options)
% [x, ~] = PwgOptimiser(xs, p, options)
% Robust Pair-Wise Geometry Optimiser
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot
%

% Configurations
switch_config;

%clf;
%plot_bundle_state(x,p);

terr=0.001;
aerr=0.0002;
delta=[terr, aerr];
gmin=0.30;
gmax=0.95;
maxitr=options.maxitr;

GATE_PROB_RESID = [0.05 0.9 0.7];

%GATE_PROB_RESID=gmin+rand(1,maxitr)*(gmax-gmin);

%GATE_PROB_RESID = rand(1, maxitr);
%GATE_PROB_RESID = ones(1,maxitr) * .5;
%GATE_PROB_RESID = min(GATE_PROB_RESID, .95); % max. prob. limit
%GATE_PROB_RESID = max(GATE_PROB_RESID, .05); % min. prob. limit
%GATE_PROB_RESID(1) = max(GATE_PROB_RESID(1),.5); % min. first gate prob. limit
if isfield(options,'save')
	save(strcat(options.save,'/gate.mat'),'-v7.3','GATE_PROB_RESID');
	load(strcat(options.save,'/gate.mat'));
end

x = xs;

i = 0;
%for i = 1 : length(GATE_PROB_RESID)
while 1
    
    i = i+1;
    x0 = x; sw0 = sw;
    
    % termination criterion
    if i>length(GATE_PROB_RESID)
        disp('Max number of iterations reached.');
        break
    end
    
    disp('    ');
    disp(['#', num2str(i), ': Minimum gate = ', num2str(GATE_PROB_RESID(i))]);
    options = constraint_graph_optimise_set_options(options, ...
        'gateresid', chi_square_bound(GATE_PROB_RESID(i), 2));
    [y, Y, sw, x] = optimise_constraint_image_inverse_depth(Ct, C, sw, x, options, ncams);
    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x(7:9)', x(10:12)'*180/pi);

    % termination criterion
    delta(1) = norm(x(7:9)-x0(7:9));
    delta(2) = norm(x(10:12)-x0(10:12));
    disp(['Update norm = ', num2str(delta)]);
    if i>1 && delta(1)<terr && delta(2)<aerr % && sum(sw)>sum(sw0)
        %if abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))>0.05
        if abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))>0.05 || ... % the case of sufficient gate change
                (abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))<=0.05&&sum(sw0)~=sum(sw)) || ... % the case of insufficient gate change
                (sum(sw&sw0)==sum(sw) && sum(sw0)==sum(sw)) || ... % the case of no switching change
                (norm(options.xest(1:3)-x(7:9))<0.0005 && norm(options.xest(4:6)-x(7:9))<0.00035) % the case of temporal stereo consistency
            disp(['Update norm = ', num2str(delta)]);
            options.xest = x(7:12);
            break
        end
    end
    
end

[x, P] = recover_moments(y, Y);
fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x(7:9)', x(10:12)'*180/pi);
