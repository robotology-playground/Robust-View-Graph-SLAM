function [x, P, C, sw, GATE_PROB_RESID] = PwgOptimiser(C, Ct, xs, sw, ncams, options)
%[x, ~] = PwgOptimiser(xs, p, options)
% Robust Pair-Wise Geometry Optimiser
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

% Configurations
config_rswitch;
if isfield(options,'save')
	save(strcat(options.save,'/gate.mat'),'-v7.3','GATE_PROB_RESID');
	%load(strcat(options.save,'/gate.mat'));
end

%clf;
%plot_bundle_state(x,p);

x = xs;
i = 0;
while 1 % applies an optimisation exit criterion
%for i = 1 : length(GATE_PROB_RESID) % uses a given number of residuals switching iterations

    i = i+1;
    x0 = x; sw0 = sw;
    
    % termination criterion
    if i>length(GATE_PROB_RESID)
        disp('Max number of iterations reached.');
        break
    end
    
    disp('    ');
    disp(['#', num2str(i), ': Minimum gate = ', num2str(GATE_PROB_RESID(i))]);
    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x(7:9)', x(10:12)'*180/pi);
    options = set_params(options, 'gateresid', chi_square_bound(GATE_PROB_RESID(i), 2));
    [y, Y, sw, x] = optimise_constraint_image_inverse_depth(Ct, C, sw, x, options, ncams);
    fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x(7:9)', x(10:12)'*180/pi);

    % termination criterion
    delta(1) = norm(x(7:9)-x0(7:9));
    delta(2) = norm(x(10:12)-x0(10:12));
    disp(['Update norm = ', num2str(delta)]);
    if i>1 && delta(1)<terr && delta(2)<aerr % && sum(sw)>sum(sw0)
        %if abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))>0.05
        if abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))>0.05 || ... % the case of sufficient gate change
                (abs(GATE_PROB_RESID(i)-GATE_PROB_RESID(i-1))<=0.05&&sum(sw0)~=sum(sw)) || ... % the case of insufficient gate change but switching change
                (sum(sw&sw0)==sum(sw) && sum(sw0)==sum(sw)) || ... % the case of no switching change
                (norm(options.xest(1:3)-x(7:9))<0.0005 && norm(options.xest(4:6)-x(7:9))<0.00035) % the case of temporal stereo consistency
            disp(['Update norm = ', num2str(delta)]);
            options.xest = x(7:12);
            break
        end
    end
    
end

[x, P] = recover_moments(y, Y);
%P = inv(Y);
%x = Y\y;

fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x(7:9)', x(10:12)'*180/pi);
