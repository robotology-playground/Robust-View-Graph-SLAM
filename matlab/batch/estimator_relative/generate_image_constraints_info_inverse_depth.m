function C = generate_image_constraints_info_inverse_depth(C, xs, ncams)
%C = generate_image_constraints_info_inverse_depth(C, xs, ncams)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

k = 0;

for i = 1:length(C)
    
    i1 = (C(i).cam-1)*6 + (1:6); % pose index
    i2 = C(i).kpt + ncams*6; % point (also inverse depth) index
    idx = [i1, i2];
    
    C(i).y = zeros(7,1);
    C(i).Y = zeros(7,7);
        
    % Numerical checks (remove zero inverse depth)
    if xs(i2,1) == 0
        k = k+1;
        continue;
    end
    
    % Compute Jacobian and {y,Y} factors for C(i)
    Hs = observation_model_jacobian_inverse_depth(xs(idx,1), C(i).p1);
    [y, Y] = canonical_update_linearised(zeros(7,1), zeros(7), ...
        @observation_model_inverse_depth, [], C(i).z, C(i).R, xs(idx,1), ...
        Hs, 1:7, [], C(i).p1);
    
    % Numerical checks (remove large information)
    if max(diag(Y)) > 1e+90
        k = k+1;
        continue;
    end
    
    C(i).y = y;
    C(i).Y = Y;
    
end

if k > 0
    if k > 1
        disp([num2str(k), ' constraints were suppressed']);
    else
        disp([num2str(k), ' constraint was suppressed']);
    end
end