% lmatch_resolve  Resolves a set of tentative matches to (sub)optimal set
% of consistent final matches.
%
% [MM,ll] = lmatch_resolve(M,l,P,I,D,opt), where the parameters mean:
%
% M (on input) ... input tentative matchces
% l, P, D ... see 'lmatch_generate'
% opt ... options, see lmatch_options
% MM ... resulting consistent matches (MM.li are indices to ll; matches are sorted by MM.C)
% ll ... resulting image line segments (ll is different from l if opt.Merging==1)

function [M,l] = lmatch_resolve(M,l,P,I,D,opt)

% Print state
fprintf('RESOLVING INCONSISTENCIES IN TENTATIVE MATCHES\n');
tic
%fprintf(' Number of views             [views                   ]: %i\n',length(P));
% s = {'OFF','ON'};
% fprintf(' Ordering constraint         [OFF|ON                  ]: %s\n',s{double(opt.Ordering)+1});
% fprintf(' Merging fragmented segments [OFF|ON                  ]: %s\n',s{double(opt.Merging)+1});
% if opt.Merging
% fprintf(' Merging residual            [pixels                  ]: %i\n',opt.MergeResidual);
% s = {'projective','affine','metric'};
% if opt.Calibration<1 | opt.Calibration>3, error('Wrong option: Calibration'); end
% fprintf(' Calibration                 [%s|%s|%s]: %s\n',s{:},s{opt.Calibration});
% fprintf(' Max. reproj. residuals      [pixels_lin pixels_nonlin]: [%g %g]\n',opt.ReprojResidual);
% fprintf(' Correl. window half size    [pixels_perp pixels_paral]: [%g %g]\n',opt.CorrParams(1),opt.CorrParams(2));
% fprintf(' Correl. window distance     [pixels                  ]: %g\n',opt.CorrParams(3));
% s = {'center','center_and_sides'};
% fprintf(' Correl. window placement    [%s|%s ]: %s\n',s{:},s{opt.CorrParams(4)+1});
% end
fprintf(' Histogram of matches [#views_per_match:#matches]: '); lmatch_print_stat(M.li);

for k = 1:length(P)
  for n = 1:length(l{k}) % precompute straight lines
    l{k}(n).l = cros(hom([l{k}(n).u l{k}(n).v]))';
  end
  if ~isfield(l{k},'s') % Compute covariance matrices of line segments
    for n = 1:length(l{k})
      uv = hom([l{k}(n).u l{k}(n).v]);
      l{k}(n).s = vgg_vech(uv*uv');
    end
  end
end

% resolve matches' inconsistency by stable matching algorithm
fprintf('\n Resolving ');
i = []; % list of final matches
ii = 1:size(M.li,2); % list of matches compatible with i still to process
vgg_hash([length(ii) 0]);
while ~isempty(ii)
  i0 = ii(argmax(M.C(ii)));  % find best match among ii
  i = [i i0];
  ii = setminus(ii,i0);
  ii = setminus( ii, inhibition_set(i0,ii,M,l,P,I,D,opt) );
  vgg_hash(length(ii));
end

M.li = M.li(:,i);
M.C = M.C(i);
M.X = M.X(:,i);
M.Y = M.Y(:,i);

fprintf(' '); lmatch_print_stat(M.li);

% remove precomputed straight lines
for k = 1:length(l)
  l{k} = rmfield(l{k},'l');
end

% merge subsets of matches having a common segment
if opt.Merging % the if-condition is not necessary, only speeds things up
  fprintf('\n Merging:');
  n = 0; % number of line segments in matches before merging
  for k = 1:length(l)
    n = n + length(unique(M.li(k,M.li(k,:)>0)));
  end
    
  [M,l] = lmatch_merge(M,l,P,cellsizes(I));
  fprintf(' '); lmatch_print_stat(M.li);
  
  m = 0; % number of line segments in matches after merging
  for k = 1:length(l)
    m = m + length(unique(M.li(k,M.li(k,:)>0)));
  end
  fprintf('\n          %i line segments merged to %i\n',n,m);
  
  % sort matches according to M.C
  [dummy,i] = sort(-M.C);
  M.li = M.li(:,i);
  M.C = M.C(i);
  M.X = M.X(:,i);
  M.Y = M.Y(:,i);
end

fprintf('\n Elapsed time: %g sec\n\n',toc);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes set of matches among ii not compatible with the match i0.
% M, l, P, I, D, opt ... as usual
% i ... matches from ii incompatible with i0 (a subset of ii)
% Index vectors i0, ii, i are all indices to M.li.
function i = inhibition_set(i0,ii,M,l,P,I,D,opt)

k0 = find(M.li(:,i0)>0);

% i := matches among ii having a common segment (in at least one view) with match i0
i = ii(any( M.li(k0,i0*ones(1,length(ii)))==M.li(k0,ii), 1 ));

% add matches violating ordering constraint with match i0
if opt.Ordering>0
  io = ordering_constraint(i0,ii,M,l,P,opt.Ordering);
  i = unique([i io]);
end

% subtract matches that have a common segment but are mergeable with li0
if opt.Merging
  if ~isempty(i)
    for i1 = i
      if lmatch_mergeable([i0 i1],M,l,P,I,D,opt)
        i(i==i1) = [];
      end
    end
  end
end

return