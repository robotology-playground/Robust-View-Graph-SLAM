%vgg_linesegs_from_edgestrips  Breaks edgel strips and approximates them by line segments.
%   [u,v] = vgg_linesegs_from_edgestrips(e) Given a set of connected edgels detected in
%   an image (eg output of vgg_xcv_segment('canny',...)), fits line segments to them.
%   e is cell vector of double(2,?), set of edgel strips; e{i}(1:2,n) is n-th edgel
%   of i-th strip. u,v double(2,?) are end points of detected line segments.
%
%   [u,v] = vgg_linesegs_from_edgestrips(e,name1,value1,...) takes optional breaking/fitting
%   parameters. E.g., vgg_linegs_from_edgestrips(e,'ChisqSigma',1,'WormLength',10).
%   The options can be :-
%      'WormLength' ... length of the worm sliding along an edge strip (default: 5 pixels)
%      'WormThreshold' ... residual threshold when the worm breaks the strip (default: .5 pixels)
%      'ChisqSigma' ... std of chi-square statistic based decision whether an edgel set
%         belongs to a single line segment (default: .45 pixels)
%      'ChisqThreshold' ... threshold of chi-square likelihood of this (default: .1)
%
%    [s,C] = vgg_linesegs_from_edgestrips(e) also returns covariance matrices of the segments.
%    C is double(6,:), C(:,n) are 6 elements of cov. matrix of segment n.
%    Get the full matrix as vgg_vech(C(:,n)).

% werner@robots.ox.ac.uk, 2001

function [u,v,C] = vgg_linesegs_from_edgestrips(e,opt)

% set default fitting parameters - see helps for the local functions below
if exist('opt')~=1
  opt = [];
end
opt = setfield_cond(opt,'WormLength',5);
opt = setfield_cond(opt,'WormThreshold',.5);
opt = setfield_cond(opt,'ChisqSigma',.45);
opt = setfield_cond(opt,'ChisqThreshold',.1);

u = [];
v = [];
C = [];
for n = 1:length(e)
  e{n} = e{n}(1:2,:);
  [ij,ijn] = chain2segments(e{n},opt.WormLength,.5);
  c = [];  
  for j = 1:size(ij,2)
    c{j} = e{n}(:,ij(1,j):ij(2,j));
  end
  c = segments2lines(c,opt.ChisqSigma,opt.ChisqThreshold);
  for j = 1:length(c)
    [u(:,end+1),v(:,end+1),Cj] = vgg_lineseg_from_x(c{j});
    C(:,end+1) = vgg_vech(Cj);
  end
end

return


%%%%%%%%%%%%%% local functions %%%%%%%%%%%%%%%%

% [ij,ijn] = chain2segments(x,W,th)  Segments chains of edgels for polygonal approximation.
% 
% x ... size (2,N), edgel coordinates
% W ... scalar, radius of the kernel sliding on the edge chain (recommended: 5)
% th ... scalar, segmentation threshold (recommended: .5)
% ij ... size (2,M), indices to segments. For each n, x(:,ij(1,n):ij(2,n)) is am edgel set of a hypothesized line segment.
% ijn ... indices to the rejected segments; same format as ij
function [ij,ijn] = chain2segments(x,W,th)

N = size(x,2);
if N < 2
  ij = [];
  return
end
px = vgg_get_homg(x);  

t = logical(zeros(1,N));

% sliding segment
r = nan*ones(1,N);
for i = 1+W : N-W
  l = px(:,i+W)'*vgg_contreps(px(:,i-W));
  d = l/sqrt(sum(l(1:2).^2))*px(:,i-W:i+W);
  r(i) = max(abs(d));
  if r(i) < th
    t(i-W+1:i+W-1) = 1;
  end
end

ij = findcontig(t);
if ~isempty(ij)
  ij = ij(:,ij(2,:)-ij(1,:)>=3); 
end
if nargout > 1
  ijn = findcontig(~t);
end
return


% sigma ... sigma of chi-square based decision whether an edgel set belongs to a line segments
% Qth ... threshold of this chi-square probability
function e = segments2lines(e,sigma,Qth)

E = []; c = 1;
for n = 1:length(e)
  if size(e{n},2) < 4
    continue
  end
  
  % try to fit straight line
  [l,r,Q] = hplanefit(e{n},sigma);

  if Q > Qth % it fits, let's put it in the output list
    E{c} = e{n};
    c = c + 1;
  
  else % it doesn't fit; let's try to break it in two
    l = vgg_wedge(vgg_get_homg(e{n}(:,[1 end]))); % line joining the segment's end points
    [dummy,i] = max(abs(l*vgg_get_homg(e{n}))); % edgel in which the residual is greatest
    y{1} = e{n}(:,1:i-1); % divide the chain in this point
    y{2} = e{n}(:,i+1:end);
    for j = 1:2 % try to fit either piece to a line
      if size(y{j},2) < 4
	continue
      end
      [l,r,Q] = hplanefit(y{j},sigma);
      if Q > Qth
        E{c} = y{j};
        c = c + 1;
      end
    end
  end
end
e = E;

% merge neighboring line segments if possible
cont = 1;
while cont
  cont = 0;
  for n = 1:length(e)-1
    [l,r,Q] = hplanefit([e{n} e{n+1}],sigma);
    if Q > Qth
      e{n} = [e{n} e{n+1}];
      e = {e{1:n} e{n+2:end}};
      cont = 1;
      break
    end
  end
end

return



function [p,r,Q] = hplanefit(x,sigma)
% [p,r,Q] = hplanefit(x [,sigma])  Fits hyperplane to points by minimizing SSD of Euclidean distances of the points from the hyperplane.
%
% x ... size (D,N), points in inhomog coords; D ... dimension, N ... number of points
% p ... size (1,D+1), hyperplane in homog coords, its normal vector has unit norm
% Q ... scalar, probability that x fits the line model, assuming noise with sigma in x (chi-square fit). For a good fit, it should be Q>.1
% r ... sum of squared residuals

x(:,all(isnan(x))) = [];
if prod(size(x)) == 0
  p = [nan nan nan nan]; 
  r = nan;
  Q = nan;
  return
end

[D,N] = size(x);
c = sum(x,2)/N;
x = x - c*ones(1,N);

[u,s,v] = svd(x*x');
a = [zeros(1,D-1) 1]*u';
p = [a -a*c];
r = s(D,D); % sum of squared residuals

if nargout > 2
  chisq = r/sigma^2;
  Q = 1 - gammainc(chisq/2,(N-D)/2); % = complementary chi-square distribution Q(chisq,N-D)
end

return


% c = findcontig(x)  Finds contiguous pieces of 1's in a logical vector x
%
% c ... size (2,?), c = [beg1 end1; beg2 end2; ...]', beg, end  are start and end indices of the found pieces
function c = findcontig(x)

x = x(:)';

d = x(1:end-1)-x(2:end);
d(end+1) = x(end);
d = [-x(1) d];

c = [find(d==-1); find(d==1)-1];

return


% Conditional setfield: sets only if the field is not present.
function s = setfield_cond(s,f,x)
if ~isfield(s,f)
  s = setfield(s,f,x);
end
return