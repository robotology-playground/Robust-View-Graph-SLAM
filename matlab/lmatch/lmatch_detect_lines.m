%[u,v] = lmatch_detect_lines(I,minLength) Detects line segments in image I.
% minLength ... minimal length of lines detected
% Works only in Windows.

function [u,v] = lmatch_detect_lines(I,minLength)

e = vgg_xcv_segment(I,'canny_edges','-sigma 1 -edge_min 20');
for n = 1:length(e)
  e{n} = e{n}([2 1],:);
end
[u,v] = vgg_linesegs_from_edgestrips(e);

% keep lines longer than a threshold
i = u-v; i = find(sum(i.*i,1) >= minLength^2);
u = u(:,i);
v = v(:,i);

return
