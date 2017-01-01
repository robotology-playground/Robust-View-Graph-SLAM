function C = constraint_loop_statistics(C)

% By Tim Bailey
% changes: by Tariq, modified to relative to relative_w.
% changes: by Tariq, modified to global to global_w.
% changes: by Tariq, from 3D case [0;0;0] to 6D [0;0;0;0;0;0].
% changes: by Tariq, added a try/catch statements to take care of
% non-square graphs.

for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = 0;
        C(i,j).no = 0;
    end
end

m = size(C,1);
n = size(C,2);
for i = 1:m
    for j = (i+1):n
        if isempty(C(i,j).z), continue, end
        for k = (j+1):n
            try
                if isempty(C(j,k).z) || isempty(C(i,k).z), continue, end
                ki = transform_to_relative_w([0;0;0;0;0;0], C(i,k).z);
                ji = transform_to_global_w(ki, C(j,k).z);
                ii = transform_to_global_w(ji, C(i,j).z);
                if sum(ii(1:3).^2) < (3.0)^2 && sum(ii(4:6).^2) < (.5*pi/180)^2
                    C(i,j).yes = C(i,j).yes + 1;
                    C(j,k).yes = C(j,k).yes + 1;
                    C(i,k).yes = C(i,k).yes + 1;
                else
                    C(i,j).no = C(i,j).no + 1;
                    C(j,k).no = C(j,k).no + 1;
                    C(i,k).no = C(i,k).no + 1;
                end
            catch
                continue
            end
        end
    end
end

%
%

% function [z, R, H] = triangle_loop(C, i, j, k)
% ki = transform_to_relative([0;0;0], C(i,k).z);
% ji = transform_to_global(ki, C(j,k).z);
% ii = transform_to_global(ji, C(i,j).z);
%
% %
%
% function z = triangle_loop_model(x)
% ki = transform_to_relative([0;0;0], C(i,k).z);
% ji = transform_to_global(ki, C(j,k).z);
% ii = transform_to_global(ji, C(i,j).z);
%
%
% function v = triangle_loop_norm(v)
% v(3,:) = pi_to_pi(v(3,:));
