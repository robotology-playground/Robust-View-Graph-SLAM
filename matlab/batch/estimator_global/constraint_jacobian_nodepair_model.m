function H=constraint_jacobian_nodepair_model(xs,idx)
% Complete the Jacobian for relative transform constraint of x2 wrt x1.
% Most of the Jacobian is based only on base frame x1, and was precomputed
% in precompute_pernode_info. This function fills in H(1:2,3).
i1=get_state_index(idx(1));
i2=get_state_index(idx(2));
H=relpose_jacobian_6DoF(xs(i1),xs(i2));
