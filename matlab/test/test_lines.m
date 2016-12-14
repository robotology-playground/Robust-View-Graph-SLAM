clear all; close all; clc;

folder = '/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab';
image = [folder, '/lmatch/examples/kampa/001.jpg'];
addpath([folder, '/lmatch']);
addpath([folder, '/lmatch/privatee']);
root = '/home/tabuhashim';
koroibot = '/Documents/MATLAB/koroibot';
addpath([root,koroibot,'/mexopencv']);

I = imread(image);
% I = rgb2gray(I);
I = uint8(sum(double(I),3)/3);

% imshow(I); drawnow;

% extract edges
edgeim = edge(I,'canny_old', .1, sqrt(2));
% link edges
[edgelist, labelededgeim] = edgelink(edgeim, 20);
drawedgelist(edgelist, size(I), 1, 'rand', 2); axis off; drawnow;
for n = 1:length(edgelist)
  e{n} = edgelist{n}';
  e{n} = e{n}([2 1],:);
end

% % extract edges
% edges = cv.Canny(I, [50 150], ...
%     'L2Gradient', 0, ...
%     'ApertureSize', 3); % edges image
% imshow(edges);

% % fit line segments to the edgelists
% tol = .1; % standard deviation from the edge in pixels
% seglist = lineseg(edgelist, tol);
% drawedgelist(seglist, size(I), 1, 'rand', 3); axis off; drawnow;

% fit line segments to the edgelists
[u,v] = vgg_linesegs_from_edgestrips(e);

% % fit line segments to the edgelists
% lines = cv.HoughLinesP(edges, ... 
%     'Rho', 1, ...
%     'Theta', pi/180, ...
%     'MaxLineGap', 0, ...
%     'Threshold', 200, ...
%     'MinLineLength', 20);
% % 'Rho', 10, ...
% %     'Theta', pi/180, ...
% %     'Threshold', 80, ...
% %     'MinLineLength', 20, ...
% %     'MaxLineGap', 10);
% for ii=1:size(lines,2);
%     u(1:2,ii) = lines{ii}([1 2])';
%     v(1:2,ii) = lines{ii}([3 4])';
% end

% keep lines longer than a threshold
i = u-v; 
i = find(sum(i.*i,1) >= 20^2);
u = u(:,i);
v = v(:,i);

% plot
s = [u; v]';
figure; imshow(I); hold on;
%plot(u(1,:), u(2,:), '.');
%plot(v(1,:), v(2,:), '.');
    
for iii=1:size(s,1);
    line([s(iii,1) s(iii,3)],[s(iii,2) s(iii,4)],'linewidth',2);
end


axis tight; drawnow;

