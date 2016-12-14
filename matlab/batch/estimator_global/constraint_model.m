function xr=constraint_model(x)
%xr=transform_to_relative(x(4:6),x(1:3)); % Tim's 3DoF
xr=transform_to_relative_w(x(7:12),x(1:6));
%xr=relpose(x(1:6),x(7:12));

