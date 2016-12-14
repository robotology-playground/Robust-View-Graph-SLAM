function plot_graph_constraints(pose, C, keep, s_1, s_2)

son = zeros(3,0); ton = son; son2 = son;
soff = son; toff = son; sof2 = son;
for i = 1:length(C)
    edge = C(i).edge;
    xs = pose(1:3, edge(1));
    
    C(i).z(1:3) = s_2*C(i).z(1:3)/norm(C(i).z(1:3)); % plot as directions
    
    xf = transform_to_global_w(C(i).z, pose(:, edge(1)));
    xf = xf(1:3);
    
    if keep(i) == 1
        son = [son xs];
        ton = [ton xf];
        son2 = [son2 pose(1:3, edge(2))];
    elseif keep(i) == 0
        soff = [soff xs];
        toff = [toff xf];
        sof2 = [sof2 pose(1:3, edge(2))];
    end
end
pon  = line_plot_conversion_(son,  ton);
pon2 = line_plot_conversion_(son2, ton);
poff = line_plot_conversion_(soff, toff);
pof2 = line_plot_conversion_(sof2, toff);

% plot3(pof2(1,:), pof2(3,:), pof2(2,:), 'b:', ... 
%      poff(1,:), poff(3,:), poff(2,:), 'm', ...
%      pon(1,:),  pon(3,:),  pon(2,:),  'g', ...
%      pon2(1,:), pon2(2,:), pon2(3,:), 'r');

plot3(pon(1,:), pon(3,:), pon(2,:), 'b', ...
     poff(1,:),poff(3,:),poff(2,:), 'r');

% plot(pon(1,:), pon(3,:), 'g', ...
%      poff(1,:),  poff(3,:),  'm');
 
hold on, axis equal, grid on

% plot camera pose frames
frame = s_1*[-0.1 1 nan,    0 0 nan, nan    0 0;
           0 0 nan, -0.2 1 nan, nan    0 0;
           0 0 nan,    0 0 nan, nan -0.2 1];
for i = 1:size(pose,2)
   fg = transform_to_global_w(frame, pose(:,i));
   plot3(fg(1,:),fg(3,:),fg(2,:), 'k')
   %plot(fg(1,:),fg(3,:), 'k')
end

% plot3(pose(1,1:2:end),pose(3,1:2:end),pose(2,1:2:end),'g');
% plot3(pose(1,2:2:end),pose(3,2:2:end),pose(2,2:2:end),'m');
plot3(pose(1,:),pose(3,:),pose(2,:),'g');