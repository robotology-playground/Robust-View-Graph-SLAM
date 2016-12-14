% lmatch_generate  Generate tentative matches based on a single image basepair.
%
% M = lmatch_generate(l,P,I,D,kk,opt), where the parameters mean:
%
% l ... line segments in images, in the following form:
%   - l{k}(n).u, l{k}(n).v ... 2-by-1 vectors, end points of n-th segment in k-th image
%   - l{k}(n).s ... 6-by-1 vector (optional), vech of inverted covariance matrix of the segment, normalised such that s(6)=2.
% P ... camera matrices, P{k} is 3-by-4 matrix of camera in k-th image. P-matrices must have correct overall scales,
%   representing an oriented projective reconstruction, i.e., there exist homogeneous 3D points X and scalars s>0
%   such that  s*[u;1]=P*X for all (non-homogeneous) endpoints u.
% I ... images, I{k} is k-th image, 2-dimensional uint8 array (only grayscale; color images are not supported)
% D ... array #views-by-#views, D(k,l) is a "distance" of views k and l.
%   It is a heuristic, saying that "nearer" views have more matches.
% kk ... 2-vector, base view pair
% opt ... option structure, see lmatch_options
% M ... set of tentative matches:
%   - M.li = line matches (indexes to l; 0 means a missing segment)
%   - M.C  = scores
%   - M.X, M.Y  = endpoints of 3D line segments (use 'lmatch_reconstruct' to make them more accurate and/or constrained)

function M = lmatch_generate(M,l,P,I,D,kk,opt)

if isempty(M) % initialize list M of tentative matches
  M.li = int32(zeros(length(P),0)); % Meaning of M.li(k,i):
                                    %  >0 = line segment in view k of match i, index to l{k}
                                    %   0 = match i has missing line segment in view k
                                    %  -1 = view k of match i has not yet been searched
  M.C = [];
  M.X = [];
  M.Y = [];
end

% print out initial state and options
fprintf('GENERATING TENTATIVE MATCHES BASED ON VIEW PAIR [%i %i]\n',kk);
tic
% for k = 1:length(P), fprintf('%i ',length(l{k})); end
% fprintf(' Number of views           [views                   ]: %i\n',length(P));
% fprintf(' Numbers of line segments  [line segments           ]: ');
% for k = 1:length(P), fprintf('%i ',length(l{k})); end
% fprintf \n
% fprintf(' Displacement search range [pixels                  ]: %g\n',opt.DispRange);
% fprintf(' Maximal view distance     [views                   ]: %i\n',opt.MaxViewDistance);
% s = {'projective','affine','metric'};
% if opt.Calibration<1 | opt.Calibration>3, error('Wrong option: Calibration'); end
% fprintf(' Calibration               [%s|%s|%s]: %s\n',s{:},s{opt.Calibration});
% fprintf(' Max. reproj. residuals    [pixels_lin pixels_nonlin]: [%g %g]\n',opt.ReprojResidual);
% fprintf(' Correl. window half size  [pixels_perp pixels_paral]: [%g %g]\n',opt.CorrParams(1),opt.CorrParams(2));
% fprintf(' Correl. window distance   [pixels                  ]: %g\n',opt.CorrParams(3));
% s = {'center','center_and_sides'};
% fprintf(' Correl. window placement  [%s|%s ]: %s',s{:},s{opt.CorrParams(4)+1});

for k = 1:length(P)
  for n = 1:length(l{k})
    l{k}(n).l = norml(cros(hom([l{k}(n).u l{k}(n).v]))'); % precompute straight lines
  end
  if ~isfield(l{k},'s') % precompute covariance matrices of line segments from endpoints, if not present
    for n = 1:length(l{k})
      uv = hom([l{k}(n).u l{k}(n).v]);
      l{k}(n).s = vgg_vech(uv*uv');
    end
  end
end

% Match the base view pair kk
fprintf(' Histogram of matches [#views_per_match:#matches]: '); lmatch_print_stat(M.li);
fprintf('\n Matching base view pair ');
if D(kk(1),kk(2)) > opt.MaxViewDistance
  error('Distance of base pair is greater than maximal view distance');
end
M = lmatch_basepair(M,kk,l,P,I,opt);
M = remove_overlaps(M);
if opt.Calibration > 1
  M = remove_beyond_infty(M,P);
end
fprintf(' '); lmatch_print_stat(M.li); fprintf \n;

% Find support in other than base views
while 1
  ui = find(any(M.li==-1,1)); % matches having yet unsearched views
  if isempty(ui) % If all views in all matches have been searched, stop
    break
  end
  fprintf(' Matching other views    ');
  vgg_hash([1 length(ui)]);
  for uii = 1:length(ui)
    i = ui(uii); % index of current match
    K0 = find(M.li(:,i)==-1); % views not yet searched
    K1 = find(M.li(:,i)>0); % views having line segments
    
    % Compute [k0 k1] = nearest views between sets K0 and K1.
    % (Note: This may be costly for large number of views, because D is large.
    % But for numbers of views manageable by the software, it turns out to be neglectable.)
    [minv,k0] = min(D(K0,K1),[],1);
    [minv,k1] = min(minv,[],2);
    k0 = K0(k0(k1));
    k1 = K1(k1);

    M.li(k0,i) = 0;
    if D(k0,k1) <= opt.MaxViewDistance
      M = lmatch_search_view(M,k0,k1,i,l,P,I,opt);
    end

    vgg_hash(uii);
  end  

  M = remove_overlaps(M);
  if opt.Calibration > 1
    M = remove_beyond_infty(M,P);
  end
  fprintf(' '); lmatch_print_stat(M.li); fprintf \n;
end
fprintf(' Elapsed time: %g sec\n\n',toc);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove matches that are subset of another match.
% (This happens only rarely.)
function M = remove_overlaps(M)
%fprintf('Removing overlaps .. ');
r = []; % matches to be removed
for i = 1:size(M.li,2)
  ki = find(M.li(:,i)>0);
  j = find( all(M.li(ki,:) == repmat(M.li(ki,i),[1 size(M.li,2)])) );
  if length(j)==1
    continue
  end
  for j = j
    if nnz(ki) < nnz(M.li(:,j)>0) % make it that we always remove match j
      j = i;
    end
    r = [r j];
  end
end
M.li(:,r) = [];
M.C(r) = [];
M.X(:,r) = [];
M.Y(:,r) = [];
%fprintf('%i removed\n',nnz(r));
return

% Remove matches that have either endpoint beyond the plane at infty (if Calibration>1)
function M = remove_beyond_infty(M,P)
%fprintf('Removing matches beyond infinity .. ');
r = logical([]);
for i = 1:size(M.li,2)
  r(i) = logical(0);
  for k = find(M.li(:,i)>0)'
    if nnz([0 0 1]*P{k}*hom(nhom([M.X(:,i) M.Y(:,i)])) < 0)
      r(i) = logical(1);
      break
    end
  end
end
M.li(:,r) = [];
M.C(r) = [];
M.X(:,r) = [];
M.Y(:,r) = [];
%fprintf('%i removed\n',nnz(r));
return
