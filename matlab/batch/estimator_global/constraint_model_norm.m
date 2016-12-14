function x = constraint_model_norm(x)
% Normalise the angle of the relative pose
x(5) = pi_to_pi(x(5));
