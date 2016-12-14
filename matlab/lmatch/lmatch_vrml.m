%lmatch_vrml  Stores reconstructed lines to VRML file.

function lmatch_vrml(filename,X,Y,d)

L = {};
for n = 1:size(X,2)
  L{n} = nhom([X(:,n) Y(:,n)]);
end
vgg_vrml_open(filename);
Vcolor = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 1 1 1]';
if nargin<4
  vgg_vrml_polylines(filename,L,'diffuseColor',[1 1 1]);
else
  vgg_vrml_polylines(filename,L,Vcolor(:,d+1));
end
vgg_vrml_close(filename);

return

%%%%%%%%%%%%%%%%%%%

function vgg_vrml_open(file)
if ischar(file)
  f = fopen(file,'w');
  if f == -1, error(['Cannot open VRML file ' file]); end
else
  f = file;
end
fprintf(f,'#VRML V2.0 utf8\n\n');
fprintf(f,'Group { children [\n\n');
if ischar(file)
  fclose(f);
end
return

function vgg_vrml_close(file)
if ischar(file)
  f = fopen(file,'a');
  if f == -1, error(['Cannot open VRML file ' file]); end
else
  f = file;
end
fprintf(f,'\n] }');
if ischar(file)
  fclose(f);
end
return

function vgg_vrml_polylines(file,L,varargin)
if ischar(file)
  f = fopen(file,'a');
  if f == -1, error(['Cannot open VRML file ' file]); end
else
  f = file;
end
if ~iscell(L)
  L = {L};
end
fprintf(f,'Shape {\n');
fprintf(f,' geometry IndexedLineSet {\n');
fprintf(f,'  coord Coordinate {\n');
fprintf(f,'   point [\n');
fprintf(f,'    %g %g %g,\n',[L{:}]);
fprintf(f,'   ]\n');
fprintf(f,'  }\n');
fprintf(f,'  coordIndex [\n');
c = 0;
for n = 1:length(L)
  nn = size(L{n},2);
  fprintf(f,'   ');
  fprintf(f,'%g,',c:c+nn-1);
  fprintf(f,'-1,\n');
  c = c + nn;
end
fprintf(f,'  ]\n');
if nargin > 2
  if isa(varargin{1},'double')
    fprintf(f,'  color Color {\n');
    fprintf(f,'   color [\n');
    fprintf(f,'    %g %g %g,\n',varargin{1});
    fprintf(f,'   ]\n');
    fprintf(f,'  }\n');
    fprintf(f,'  colorPerVertex FALSE\n');
  else
    fprintf(f,' }\n');
    fprintf(f,' appearance Appearance {\n');
    fprintf(f,'  material');
    vgg_vrml_material(f,varargin{:});
  end
end
fprintf(f,' }\n');
fprintf(f,'}\n\n');
if ischar(file)
  fclose(f);
end
return

function vgg_vrml_material(file,varargin)
if isempty(varargin)
  return
end
if ischar(file)
  f = fopen(file,'a');
  if f == -1, error(['Cannot open VRML file ' file]); end
else
  f = file;
end
allowed = {'ambientIntensity'
           'diffuseColor'
           'specularColor'
           'emissiveColor'
           'shininess'
           'transparency'};

m = struct(varargin{:});
fprintf(f,'  Material {\n');
for s = fieldnames(m)'
  s = s{:};
  if ~any(strcmp(s,allowed))
    error(['Unknown material property: ' s]);
  end
  fprintf(f,'   %s ',s);
  fprintf(f,'%g ',getfield(m,s));
  fprintf(f,'\n');
end
fprintf(f,'  }\n');
if ischar(file)
  fclose(f);
end
return