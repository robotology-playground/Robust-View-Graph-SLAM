function batch_lines(name, P, opt)

for k = 1:length(name)
  [dummy,name{k}] = fileparts(name{k});
  I{k} = imread([base_dir '/' name{k} '.' ext]);
  if ndims(I{k})==3
    I{k} = uint8(sum(double(I{k}),3)/3);  % convert to grayscale
  end
  imsize(:,k) = size(I{k})';
  P{k} = load([base_dir '/' name{k} '.P']);
  lines_name = [base_dir '/' name{k} '.lines'];
  if exist(lines_name) ~= 2
    fprintf('\nCOMPUTING LINE SEGMENTS IN IMAGE %s',name{k});
    [u,v] = lmatch_detect_lines(I{k},20);
    fprintf(' ... %i lines detected\n',size(u,2));
    s = [u; v]'; save('-ascii',lines_name,'s');
  end
  lk = load(lines_name)';
  for n = size(lk,2):-1:1
    l{k}(n).u = lk([1 2],n);
    l{k}(n).v = lk([3 4],n);
  end
  l{k} = lmatch_lineseg_orient(l{k},I{k}); % orient image line segments
end

% Compute table D of pairwise view distances 
% Note: The method below works only for metric calibration. Choose D based on knowledge specific to the view set.
for k = 1:length(P)
  C(:,k) = null(P{k}); % homog. camera center
end
for k1 = 1:length(P)
  for k2 = 1:length(P)
    D(k1,k2) = norm(C(1:3,k1)/C(4,k1) - C(1:3,k2)/C(4,k2));
  end
end

opt = lmatch_options;
opt.Calibration = 2;
opt.DispRange = Inf; % change according to view set !!!
opt.Ordering = 1; % change to 2 or 0 to speed up resolving
opt.Merging = 1; % change to 0 to speed up resolving

if exist([base_dir '/M.mat']) == 2
  load([base_dir '/M.mat']); % load tentative matches
else
  % Generate and save tentative matches
  M = [];
  for kk = lmatch_find_basepairs(D)
    M = lmatch_generate(M,l,P,I,D,kk,opt); % change DispRange for other examples
  end
  save([base_dir '/M.mat'],'M');
end

% Discard tentative matches with less then x views and resolve
i = sum(M.li>0)>=4; % change according to view set (eg, use 4 for 'examples/kampa') !!!
MM.li=M.li(:,i);
MM.C=M.C(i);
MM.X=M.X(:,i);
MM.Y=M.Y(:,i);
[MM,ll] = lmatch_resolve(MM,l,P,I,D,opt);

% Reconstruct
%V = eye(4,3);
%V = {V(:,1),V(:,2),V(:,3),V(:,[2 3]),V(:,[3 1]),V(:,[1 2]),[]}; % use for 'merton'
%V = {V(:,1),V(:,[2 3]),[]}; % use for 'merton'
V = {[]}; % use for all other examples
[MM,d] = lmatch_reconstruct(MM,ll,P,imsize,V,3); % can be omitted for constrained (ie, with V ~= {[]}) reconstruction

% visualise reconstructed lines
%lmatch_vrml([base_dir '/lines3d.wrl'],MM.X,MM.Y); %,d);
%eval(['! vrmlview ' base_dir '/lines3d.wrl &']);

figure; 
plot3([MM.X(1,:); MM.Y(1,:)], ...
      [MM.X(2,:); MM.Y(2,:)], ...
      [MM.X(3,:); MM.Y(3,:)],'r');
  axis equal;

return

% see matches
a = images(I);
for k = 1:length(I)
  axes(a(k));
  for n = 1:size(MM.li,2)
    i = MM.li(k,n);
    if i>0
      lki = ll{k}(i);
      plot([lki.u(2) lki.v(2)],[lki.u(1) lki.v(1)],'r');
      text(lki.u(2),lki.u(1),num2str(n),'color','y','fontsize',9,'clipping','on');
    end
  end
end