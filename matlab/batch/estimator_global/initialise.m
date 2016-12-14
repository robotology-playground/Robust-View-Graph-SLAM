function [C, Ct, sw, xs] = initialise(C)
% Configurations
config_rswitch_v2;
% Generate initial path (xs) and trusted constraints (sw==1)
if USE_MAX_SPANNING_TREE
    disp('Skeleton: Max weight spanning tree')
    [xs, sw] = pose_generate_spanning_tree(C);    
else
    disp('Skeleton: Sequential path')
    xs = pose_generate_sequential(C);
    sw = diff(vertcat(C.edge), 1, 2) <= NUM_INITIAL_PATH;
end
xs = xs(:);

% If using sequential path, find and switch off any bad constraints
if ~USE_MAX_SPANNING_TREE&&REFINE_INITIAL_PATH
    disp(['Refining initial path, with thickness ' num2str(NUM_INITIAL_PATH)])
    on = find(sw==1);
    Cinfo = generate_constraint_info_pose(C(on), xs);
    options = constraint_graph_optimise_set_options;
    options = constraint_graph_optimise_set_options(options, 'verbose', 2);
    N = length(xs)/6;
    [y, Y] = initialise_info_matrix(Cinfo, N);
    [~,~, swon, xs] = constraint_graph_optimise_with_switching(y, Y, Cinfo, ones(size(Cinfo)), xs, options);
    sw(on(swon==0)) = 0; % turn off rejected constraints
end

if TRUST_INIT_LINKS
    Ct = C(sw==1);
    C = C(sw==0);
    sw = sw(sw==0);
else 
    Ct = [];
end

if INIT_WITH_HIGH_WEIGHT_LINKS
    w = vertcat(C.w);
    if INIT_WITH_HIGH_WEIGHT_LINKS == -1
        INIT_WITH_HIGH_WEIGHT_LINKS = mean(w);
    end
    disp(['Turning on links with weight greater than: ' num2str(INIT_WITH_HIGH_WEIGHT_LINKS)])
    sw = sw | w>INIT_WITH_HIGH_WEIGHT_LINKS;
end

return

function [y, Y] = initialise_info_matrix(C, N)
[y, Y] = generate_joint_info_matrix(C, N);
Y(1:6,1:6) = eye(6) * 1e4;
%spy(Y);pause;
return