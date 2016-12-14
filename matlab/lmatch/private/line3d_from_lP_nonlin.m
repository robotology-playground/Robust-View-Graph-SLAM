% line3d_from_lP_nonlin  Non-linear estimation of 3D line segment from image line segments.
% As vgg_line3d_from_lP_nonlin but faster.
%
% SYNOPSIS
% L = line3d_from_lP_nonlin(u,v,P,imsize,L0 [,opt])
%
% u, v ... double(2,K), lines end points
% P ... K-cell with 3-by-4 camera matrices
% L0 ... double(4,2), initial scene line (optional). Homogeneous points L0(:,i) span 
%   the line. If omitted, linear estimation is done first.
% L ... double(4,2), estimated 3D line. Points L(:,i) span the line.

function L = line3d_from_lP_nonlin(u,v,P,imsize,L,varargin)

old_warning = warning;
warning off MATLAB:nearlySingularMatrix

K = length(P); % number of images

% Preconditioning
if ~isempty(imsize)
  H = vgg_conditioner_from_image(imsize);
  u = nhom(H*hom(u));
  v = nhom(H*hom(v));
  for k = 1:K
    P{k} = H*P{k};
  end
  scale = H(1,1); % save the scales for evaluating objective function
else
  scale = 1;
end

if nargin<5
  L = [];
end
if isempty(L)
  for k = 1:K
    uv = hom([u(:,k) v(:,k)]);
    s{k} = uv*uv';
  end
  L = vgg_line3d_from_lP_lin(s,P);
end

% Take first and last lines as the parameterizing ones
k = [2 K];
nk = k([2 1]);
u(:,k) = u(:,nk);
v(:,k) = v(:,nk);
P(k) = P(nk);

% L is parameterized by two planes, corresponding to image lines 1, 2
%   A1 = [1 p(1:2)']*AP(:,:,1),
%   A2 = [1 p(3:4)']*AP(:,:,2)
for k = 1:2
  [dummy,dummy,V] = svd(cros(P{k}*L)');
  AP(:,:,k) = V'*P{k};
end
Ppv = line3d_Ppv(vertcat(P{:}));

% optimization
p = newton(@obj_line3d_from_lP, {Ppv,AP,cat(3,u,v)},...
           [0;0;0;0],...
           scale,...
           varargin{:});
if ~any(isnan(p)|isinf(p))
  [dummy,dummy,L] = svd([[1 p(1:2)']*AP(:,:,1)
                         [1 p(3:4)']*AP(:,:,2)]);
  L = L(:,3:4);
else
  fprintf('Divergence of line3d_nonlin\n');
end

warning(old_warning);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = newton(RES,PARAMS,a,scale, opt_iter,opt_drmsrel)

if nargin<5, opt_iter = 10; end
if nargin<6, opt_drmsrel = 1e-6; end
lambda_max = 2;
lambda_min = -8;
verbose = 0;

lambda = lambda_min;
if verbose
  e = feval(RES,a,PARAMS{:});
  ssd = sum( (e/scale).^2 );
  fprintf('initial:  lambda=%3g  rms=%.15g\n',lambda,sqrt(ssd/length(e)));
  if ssd==0, return, end
end
for i = 1:opt_iter

  % Compute actual residual and jacobian
  [e,J] = feval(RES,a,PARAMS{:});
  ssd0 = sum( (e/scale).^2 );

  JJ = J'*J;
  Je = (e'*J)';
  cnt = 0;
  while lambda<=lambda_max
    cnt = cnt + 1;
    da = -(JJ + 10^lambda*diag(diag(JJ))) \ Je;
    %da = -JJ\Je;
    ssd = sum(( feval(RES,a+da,PARAMS{:})/scale ).^2);
    if ssd==0, return, end
    drmsrel = sqrt(ssd0/ssd)-1;
    if drmsrel>0
      lambda = max(lambda_min,lambda-2);
      break
    end
    lambda = lambda + 2;
  end
  %fprintf('%i ',cnt);

  if verbose
    fprintf('iter=%i:  lambda=%3g  rms=%.15g  drmsrel=%g\n',i,lambda,sqrt(ssd0/length(e)),drmsrel);
  end
    
  if drmsrel<0, break, end
  a = a + da;
  if abs(drmsrel)<opt_drmsrel, break, end
end
%fprintf('%i',i); if abs(drmsrel)>opt_drmsrel, fprintf !; end; fprintf(' ');
if verbose
  fprintf \n
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective function
% This one is not used; instead, mex obj_line3d_from_lP.cxx is used.
function [y,J] = F(p,Ppv,AP,u)

K = size(u,2);

A1 = [1 p(1:2)']*AP(:,:,1);
A2 = [1 p(3:4)']*AP(:,:,2);
Lpv = [cros(A1(1:3)',A2(1:3)')' A2(4)*A1(1:3)-A1(4)*A2(1:3)];

y = [];
for k = 1:K
  l = Lpv*Ppv(:,(1:3)+(k-1)*3);
  l = l/sqrt(l(1)^2+l(2)^2);
  y = [ y l*[u(:,k,1);1] l*[u(:,k,2);1] ];
end
y = y';

% else, compute also jacobian
if nargout < 2
  return
end

dif = 1e-8;
J = zeros(length(y),length(p));
for i = 1:length(p)
  pdif = p;
  pdif(i) = pdif(i) + dif;
  J(:,i) = (F(pdif,Ppv,AP,u) - y)/dif;
end
return
