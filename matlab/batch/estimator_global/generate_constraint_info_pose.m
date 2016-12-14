function C = generate_constraint_info_pose(C, xs)
k = 0;
for i = 1:length(C) % parallel for loop
    edge = C(i).edge;
    ii = [get_state_index(edge(1)) get_state_index(edge(2))];
    % Compute Jacobian and {y,Y} factors for C(i)
    Hs = constraint_jacobian_nodepair_model(xs, edge);
    [C(i).y, C(i).Y] = canonical_update_linearised(zeros(12,1), zeros(12), ...
        @constraint_model, @constraint_model_norm, C(i).z, C(i).R, xs(ii), Hs, 1:12);
    % Numerical checks : very small motion
    % Numerical checks (remove large information)
    if max(diag(C(i).Y)) > 1e+9
        k = k+1;
        C(i).y = zeros(12,1);
        C(i).Y = zeros(12);
    end
end

if k > 0
    if k > 1
        disp([num2str(k), ' constraints were suppressed']);
    else
        disp([num2str(k), ' constraint was suppressed']);
    end
end