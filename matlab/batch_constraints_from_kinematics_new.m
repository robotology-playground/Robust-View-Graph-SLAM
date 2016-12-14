function [C, P, A0] = batch_constraints_from_kinematics_new(encoders, C, options)

[P, A0] = cameras_from_kinematics(encoders.eyes, encoders.neck, encoders.waist);
for k = 1:length(C);
    i = C(k).edge(1);
    j = C(k).edge(2);
    kin = [P{i};0 0 0 1]\[P{j};0 0 0 1];
    R = kin(1:3,1:3);
    t = kin(1:3,4);
    C(k).q = R2q(R);
    C(k).r = R;
    C(k).t = t';
end

% save to .mat file (this function can be called by other than the batch implementation)
if isfield(options,'save')
    save(strcat(options.save,'/fwd_kinematics'), 'C', 'P', 'A0');
end
fprintf('Done ! \n');