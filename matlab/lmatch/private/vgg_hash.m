%VGG_HASH  Prints hashes as for-loop is being executed.
%   Calling this function inside a for-loop body causes printing hashes as the loop
%   is being executed.
%
%   VGG_HASH(b) where b is 2-vector initializes the function such that the range
%      b=[minx maxx] of loop-controlling variable corresponds to printing 25 hashes.
%   VGG_HASH(b,N) uses N rather than 25 hashes.
%   VGG_HASH(b,N,s) uses string s rather than '#'.
%   VGG_HASH(x) where x is a scalar prints appropriate number of hashes.
%   c = VGG_HASH(x) does the same but returns string rather than prints it.
%
%   EXAMPLE of usage is very easy:
%   vgg_hash([2 -2])
%   for x = linspace(2,-2,123)
%     pause(.02); % loop body
%     vgg_hash(x)
%   end
%   fprintf('\n');

% (c) T.Werner, April 2002

function ss = vgg_hash(x,NN,hh)

persistent X N n h

if all(size(x)>1)
  error('x must be scalar or 2-vector');
end

switch length(x)

 case 2 % initialize

  if nargin < 2
    N = 25;
  else
    N = NN;
  end
  X = x;
  if abs(X(1)-X(2))<=eps
    X(1) = X(2) - 1;
  end

  n = 0;
  if nargin < 3
    h = '#';
  else
    h = hh;
  end

 case 1 % step

  m = floor( N*(x-X(1))/(X(2)-X(1)) );
  s = h(ones(1,m-n));
  n = m; % n = number of hashes already printed
  if nargout > 0
    ss = s;
  else
    fprintf(s);
  end
  
 otherwise
  error('x must be a scalar or 2-vector');
  
end


return