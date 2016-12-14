function z = observation_model_inverse_depth(x, p)
%z = observation_model_inverse_depth(x, p)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% range
r = 1./x(7:end)'; % range from inverse range
%r = x(7:end)'; % range

% 3d points in camera 1 frame
rim = sqrt(p(1,:).*p(1,:) + p(2,:).*p(2,:) + 1);
d = r/rim;
p(1,:) = d.*p(1,:);
p(2,:) = d.*p(2,:);
p(3,:) = d;

% image 2 measurements
p(1,:) = p(1,:) - x(1);
p(2,:) = p(2,:) - x(2);
p(3,:) = p(3,:) - x(3);
z = pflat(w2R(x(4:6))'*p(1:3,:));
z = z(1:2, :);
z = z(:);