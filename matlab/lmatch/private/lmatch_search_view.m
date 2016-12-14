% Search view k0 for presence of segment(s) compatible with match i0.
%
% k1 ... the view nearest to k0 among views in which i0 has segments

function M = lmatch_search_view(M,k0,k1,i0,l,P,I,opt)

K = length(l);

% ki := views in which the query match has segments.
% Further, we reorder ki such that ki(end)==k1.
ki = find(M.li(:,i0)>0)';
ki(ki==k1) = [];
ki(end+1) = k1;

ki0 = [ki k0]; % views in which new matches will have segments

% query
lk = {};
for k = ki
  lk{end+1} = l{k}(M.li(k,i0));
end
[n,c] = search_view({lk{:} l{k0}},P(ki0),I(ki0),[M.X(:,i0) M.Y(:,i0)],opt);
if isempty(n)
  return
end


% Update M

i_new = size(M.li,2)+[1:length(n)-1]; % make indices of new matches at the end of M.li
ii = [i0 i_new];

M.li(:,i_new) = -1;
Ci0  = M.C(i0);  %Ci0  = M.C {i0};
for i = 1:length(ii)

  M.li(ki0,ii(i)) = [M.li(ki,i0); int32(n(i))];
  M.C(ii(i)) = Ci0 - sum(log(1-c(:,i)));  %M.C {ii(i)} = [Ci0  c(:,i)'];

  % reconstruct 3d line segment
  lsi = [lk{:} l{k0}(n(i))];
  if 1 % faster
    L = line3d_from_lP_nonlin([lsi.u],[lsi.v],P(ki0),cellsizes(I(ki0)),[],5,1e-8);
  else % slower
    for k = 1:length(lsi), s{k} = vgg_vech(lsi(k).s); end
    L = vgg_line3d_from_lP_nonlin(s,P(ki0),cellsizes(I(ki0)),[],[],'niter_term',4,'loglambda_init',-15,'verbose',0);
  end
  L = lineseg3d_from_L(lsi,P(ki0),L);
  M.X(:,ii(i)) = L(:,1);
  M.Y(:,ii(i)) = L(:,2);

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Given a query match with segments in images 1:K-1, it tries to find 
% corresponding segment(s) in image K.
% View K-1 is nearest to view K (for photommetric test).
%
% l ... cell(K), l{1:K-1} are lineseg(1), l{K} is lineseg(?)
% P ... cell(K)
% I ... cell(K)
% L ... 3d segment already estimated from segments l{1:K-1}
% nK ... double(N), indices to l{K} of segments that passed all tests
% c ... double(K-1,N), scores of view pairs (1,K)..(K-1,K)
%
function [nK,c] = search_view(l,P,I,L,opt)

K = length(l);

c = [];
nK = [];

% nK := segments within displacement search range
nK = find(all( abs(vertcat(l{K}.l)*hom([l{K-1}.u l{K-1}.v]))' <= opt.DispRange ));
if isempty(nK)
  return
end

% Leave in nK only segments in the last image that are in a band given by the query segment
nK = nK( lineseg_near(l{K}(nK),nhom(P{K}*L)) );
if isempty(nK)
  return
end

% reprojection test
b = logical([]);
Lr = {}; % 3D straight lines estimated inside reprojection test (used for photommetric score)
for n = nK
  [b(end+1),Lr{end+1}] = lmatch_reprojection_test([l{1:K-1} l{K}(n)],P,cellsizes(I),opt.ReprojResidual); %,L);
end
nK = nK(b);
if isempty(nK)
  return
end
Lr = Lr(b);

% Photometric test between image K and images 1:K-1.
% The test is done between view K and its nearest view.
for n = 1:length(nK)
  lK = l{K}(nK(n));
  cn = lmatch_photoscore([l{K-1}.u l{K-1}.v],[lK.u lK.v],P{[K-1 K]},I{[K-1 K]},opt.NCCWindow,opt.Calibration>=3,Lr{n});
  cn = cn(cn>opt.NCCThreshold(1));
  if length(cn) < opt.NCCThreshold(2)/opt.NCCWindow(3)
    c(n) = 0;
  else 
    c(n) = mean(cn);
  end
end
nK = nK(c>0);
c = c(:,c>0);

return


% n = lineseg_nearest(L,l) Finds line segments among L near segment u.
%
% L ... lineseg(?)
% u ... double(2,2), query segment
%
% Segment orientations matter.
%
function n = lineseg_near(L,u)

d = normx(u(:,2)-u(:,1))';

% discard all segments that have very different orientation
n = find( d*normx([L.v]-[L.u]) > 0.75 );
if isempty(n)
  return
end

n( -d*[eye(2) -u(:,2)]*hom([L(n).u]) < 0 ) = [];
if isempty(n)
  return
end

n(  d*[eye(2) -u(:,1)]*hom([L(n).v]) < 0 ) = [];

return