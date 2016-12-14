function plot_pose_graph_w(x, i)
cols = 'krbgcmy';
persistent IDX
if isempty(IDX) || IDX == length(cols), IDX = 0; end
IDX = IDX + 1;
if nargin == 2, IDX = mod(i-1, length(cols)) + 1; end

plot3(x(1:6:end), x(3:6:end), x(2:6:end), cols(IDX))
hold on

frame = .001*[-0.1 1 nan,    0 0 nan, nan    0 0;
            0 0 nan, -0.2 1 nan, nan    0 0;
            0 0 nan,    0 0 nan, nan -0.2 1];
idx = 1:6;
while idx(1) < length(x)
    fg = transform_to_global_w(frame, x(idx));
    idx = idx + 6;
    plot3(fg(1,:),fg(3,:),fg(2,:), cols(IDX))
end
axis equal, grid on
drawnow
