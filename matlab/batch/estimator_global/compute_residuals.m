function r = compute_residuals(x, P, C, xs)

r = zeros(size(C));
for i = 1:length(C)
    edge = C(i).edge;
    ii = [get_state_index(edge(1)) get_state_index(edge(2))];
    
    % Constraint model Jacobian
    Hs = constraint_jacobian_nodepair_model(xs, edge);
    
    % Predicted constraint
    dx = x(ii) - xs(ii);
    %     %%% adde by Tariq, using rotation matrix to compute angular difference
    %     Rzst = w2R(xs(4:6))';
    %     Rz = w2R(x(4:6));
    %     Rv = Rzst*Rz;
    %     dx(4:6) = R2w(Rv)';
    %     Rzst = w2R(xs(10:12))';
    %     Rz = w2R(x(10:12));
    %     Rv = Rzst*Rz;
    %     dx(10:12) = R2w(Rv)';
    %     %%%
    dx([5 11]) = pi_to_pi(dx([5 11])); % FIXME: do we need this normalisation?
    zpred = constraint_model(xs(ii)) + Hs*dx;
    
    % Normalised residual
    v = constraint_model_norm(C(i).z - zpred);
    S = Hs*P(ii,ii)*Hs' + C(i).R;
    r(i) = v'*(S\v);
end
