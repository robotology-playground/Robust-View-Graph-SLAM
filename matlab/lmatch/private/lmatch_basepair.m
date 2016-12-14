function M = lmatch_basepair(M,kk,l,P,I,opt)

if opt.Calibration > 1
  A = [0 0 0 1]; % plane at infty
else
  A = [];
end
u2 = hom([l{kk(2)}.u]);
v2 = hom([l{kk(2)}.v]);
EB = epibeam_init(u2,v2,P{kk},A);
ll2 = norml(vertcat(l{kk(2)}.l));

Mli = M.li(kk,:);
li = int32([]);
C = [];
X = [];
Y = [];

vgg_hash([1 length(l{kk(1)})]);
for n1 = 1:length(l{kk(1)})

  i1 = (Mli(1,:)==n1);
  l1 = l{kk(1)}(n1);
  uv1 = hom([l1.u l1.v]);
  
  nn2 = epibeam_search(EB,uv1);
  nn2 = nn2(all( abs(ll2(nn2,:)*uv1)' <= opt.DispRange )); % discard segments outside disparity search range
  for n2 = nn2

    % check if not already in the list
    if any(i1 & (Mli(2,:)==n2))
      continue
    end

    l2 = l{kk(2)}(n2);
    c = lmatch_photoscore([l1.u l1.v],[l2.u l2.v],P{kk},I{kk},opt.NCCWindow,opt.Calibration>=3);
    c = c(c>opt.NCCThreshold(1));
    if length(c) < opt.NCCThreshold(2)/opt.NCCWindow(3);
      c = 0;
    else 
      c = mean(c);
    end

    if c > 0
      % reconstruct 3d line segment
      ln = [l{kk(1)}(n1) l{kk(2)}(n2)];
      Ln = vgg_line3d_from_lP_lin({vgg_vech(ln(1).s) vgg_vech(ln(2).s)},P(kk),cellsizes(I(kk)));
      Ln = lineseg3d_from_L(ln,P(kk),Ln);

      % check if any more end points are beyond plane at infinity
      if (opt.Calibration > 1) & any(any( [P{1}(3,:); P{2}(3,:)]*hom(nhom(Ln)) < 0 ))
	      continue
      end

      li = [li [n1;n2]];
      C = [C c];
      X = [X Ln(:,1)];
      Y = [Y Ln(:,2)];
    end
  
  end

  vgg_hash(n1);
end

% update M
if ~isempty(li)
  li(kk,:) = li;
  li(setminus(1:size(M.li,1),kk),:) = -1;
  M.li = [M.li int32(li)];
  M.C = [M.C -log(1-C)];
  M.X = [M.X X];
  M.Y = [M.Y Y];
end

return