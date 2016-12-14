function x = int_o1(x_dot, xprev, dt)

% int_o1 - 1st order (Euler) integration step

x = xprev + x_dot.*dt;