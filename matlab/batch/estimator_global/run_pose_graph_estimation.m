function [y, Y, sw, x, C, Ct] = run_pose_graph_estimation(C, Ct, sw, xs)
tic;

% Configurations
config_rswitch_v2;

assert(size(C,1)==1,'Need to process C via convert_C_to_CG()')
if nargin == 1
    %C = apply_nominal_R(C);
    [C, Ct, sw, xs] = initialise(C);
end

options = constraint_graph_optimise_set_options;
options = constraint_graph_optimise_set_options(options, 'verbose', VERBOSE);
options = constraint_graph_optimise_set_options(options, 'gateratio', K_RESID);
options = constraint_graph_optimise_set_options(options, 'iterations', 2);
options = constraint_graph_optimise_set_options(options, ...
    'gateinnov', chi_square_bound(GATE_PROB_INNOV, 6));

x = xs;

% pose = reshape(x, 6, length(x)/6);
% figure(1);
% subplot(3,1,1); plot(pose(1,:), 'r'); hold on;
% subplot(3,1,2); plot(pose(2,:), 'r'); hold on;
% subplot(3,1,3); plot(pose(3,:), 'r'); hold on;
% figure(2);
% subplot(3,1,1); plot(pose(4,:)*180/pi, 'r'); hold on;
% subplot(3,1,2); plot(pose(5,:)*180/pi, 'r'); hold on;
% subplot(3,1,3); plot(pose(6,:)*180/pi, 'r'); hold on;
% drawnow;

colors = hsv(length(GATE_PROB_RESID));
for i = 1:length(GATE_PROB_RESID)
    disp('    ');
    disp(['#',num2str(i)]);
    disp(['Minimum gate: ' num2str(GATE_PROB_RESID(i))]);
    options = constraint_graph_optimise_set_options(options, ...
        'gateresid', chi_square_bound(GATE_PROB_RESID(i),6));
    [y, Y, sw, x] = optimise_constraint_graph(Ct, C, sw, x, options); % relinearise final
    
%     pose = reshape(x, 6, length(x)/6);
%     figure(1);
%     subplot(3,1,1); plot(pose(1,:), 'color', rand(1,3)); hold on;
%     subplot(3,1,2); plot(pose(2,:), 'color', rand(1,3)); hold on;
%     subplot(3,1,3); plot(pose(3,:), 'color', rand(1,3)); hold on;
%     figure(2);
%     subplot(3,1,1); plot(pose(4,:)*180/pi, 'color', rand(1,3)); hold on;
%     subplot(3,1,2); plot(pose(5,:)*180/pi, 'color', rand(1,3)); hold on;
%     subplot(3,1,3); plot(pose(6,:)*180/pi, 'color', rand(1,3)); hold on;
%     drawnow;
    
end

toc
