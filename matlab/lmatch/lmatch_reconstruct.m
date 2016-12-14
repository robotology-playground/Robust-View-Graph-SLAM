% [X,Y,dir] = lmatch_reconstruct(l,li,P,imsize [,V,sigma])  Reconstructs 3D lines segments given image line segments matches.
%
% l ... lineseg(?){K}, image line segments
% li ... double(K,?), matches, indices to l
% P ... cell{K}, camera matrices
% imsize ... double(2,K), image sizes
%
% 3D lines can be optionally classified into direction classes and reconstructed with these direction constraints. 
% This increases accuracy. E.g., direction classes can be horizontal/vertical vanishing points and lines in architectural scenes.
% V ... (optional) cell(N), list of direction classes:
%   - If V{dir} is double(4,1), L passes through homog point V{dir}.
%   - If V{dir} is double(4,2), L passes through line spanned by homog points V{dir}(:,1) and V{dir}(:,2).
%   - If V{dir}==[], L is unconstrained.
% For each line, classes dir=1:N are taken one after the other. For each dir, 3D line L is reconstructed with direction constraint V{dir}.
% If this L has reprojection residual smaller than ReprojResid, it is accepted and returned along with the value of dir.
% If L cannot be reconstructed with the residual smaller then ReprojResid for any dir, dir=0 is returned.

function [M,dir] = lmatch_reconstruct(M,l,P,imsize,varargin)

fprintf('RECONSTRUCTING 3D SEGMENTS ');
tic
vgg_hash([1 size(M.li,2)]);
for n = 1:size(M.li,2)

  % form variables for the match n
  lin = M.li(:,n);
  ln = [];
  s = {};
  for k = find(lin)'
    ln = [ln l{k}(M.li(k,n))];
    s{end+1} = vgg_vech(ln(end).s);
  end
  Pn = P(lin~=0);
  imsizen = imsize(:,lin~=0);

  % estimate 3D line
  [L,dir(n)] = line3d_from_lP_classifydir(s,Pn,imsizen,[],varargin{:}); % surprisingly, if existing L is given as initial estimate results are worse
    
  % clip to line segment
  L = lineseg3d_from_L(ln,Pn,L);
  M.X(:,n) = L(:,1);
  M.Y(:,n) = L(:,2);
  
  vgg_hash(n);
end
fprintf('\n');
fprintf(' Elapsed time: %g sec\n\n',toc);

return