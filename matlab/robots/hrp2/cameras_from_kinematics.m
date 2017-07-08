function P = cameras_from_kinematics(encoders, floatingbase)

% No HRP2 kinematics anyway, this is dummy

P = cell(1,2*size(encoders,2));

if size(P,2) > 1;
    disp(' - Computing camera poses from HRP2 kinematics.');
else
	disp(' - No HRP2 kinematics.');
    return;
end
