% l = lineseg_orient(l,I)  Orients line segment(s) using image I such that the lighter part is on the rigt.
%
% l ... image line segment
% I ... double (:,:), image

function l = lmatch_lineseg_orient(l,I)

h = vgg_gauss_mask(1,1,-3:3);

for i = 1:size(l,1)
  for j = 1:size(l,2)
    if vgg_conv_lineseg(I,h,l(i,j).u,l(i,j).v) < 0
      w = l(i,j).u;
      l(i,j).u = l(i,j).v;
      l(i,j).v = w;
    end
  end
end

return