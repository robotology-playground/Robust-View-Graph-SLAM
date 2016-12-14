%VGG_XCV_SEGMENT  VXL-based corner/edge/line detection.
%   This function calls a binary executable of vxl code for computing
%   Harris corners, Canny edgel strips, and approximating these strips
%   with line segments.
%
%   s = VGG_XCV_SEGMENT(I,action) performs the given action in gray
%   or color image I (must be uint8). Action string can be :-
%    - 'corners' ... Harris corner detection. Then s is double(2,?) with
%      the corners coordinates.
%    - 'canny_edges' ... Canny edge detection + linking edgels to strips.
%      Then s{n}(4,?) is n-th strip and :-
%       - s{n}(1:2,:) are edgel coordinates,
%       - s{n}(3,:) are gradients,
%       - s{n}(4,:) are edgel orientations.
%    - 'break_edges' ... like canny_edges but the strips are additionally
%      broken in high curvature places.
%    - 'lines' ... break_edges + approximation with line segments if possible.
%      s is double(4,:).
%      See VGG_LINESEGS_FROM_EDGESTRIPS for better and slower code.
%   NOTE: All coordinates are in MATLAB's (1,1)-based coordinate frame.
%
%   s = VGG_XCV_SEGMENT(I,action,opt) allows changing default options passed
%   to the executable. E.g., opt='-sigma 3  -low 1  -high 4'.
%   Options for corner detection and defaults :-
%    -gauss_sigma float                0.7
%    -corner_count_max integer         300
%    -relative_minimum float           1e-005
%    -scale_factor float               0.04
%    -adaptive bool                    1
%   Options for edge/line detection and defaults :-
%    -sigma float                      1
%    -max_width integer                50
%    -gauss_tail float                 0.0001
%    -low float                        2
%    -high float                       12
%    -edge_min integer                 60
%    -min_length integer               10
%    -border_size integer              2
%    -border_value float               0.0
%    -scale float                      5.0
%    -follow_strategy integer          2
%    -join_flag bool                   1
%    -junction_option integer          0
%    -bk_thresh float                  0.3
%    -str_high float                   12
%    -str_edge_min integer             60
%    -str_min_length integer           60
%    -min_fit_length integer           10
%
%   See also VGG_LINESEGS_FROM_EDGESTRIPS, VGG_HARRIS, VGG_HARRIS_AFFINE.

% (c) binary: kym@robots.ox.ac.uk, May 2002
%     matlab function: werner@robots.ox.ac.uk, May 2002

function s = vgg_xcv_segment(I,action,opt)

if ~isa(I,'uint8')
    error('Image must be uint8');
end

inname  = 'vgg_xcv_in.png';
outname = 'vgg_xcv_out';

% Save input image
imwrite(I,inname,'png');

% Form options
if nargin < 3
    opt = '';
end
switch action
    case 'corners'
        opt = [' -gauss_sigma        0.7',...
            ' -corner_count_max   300',...
            ' -relative_minimum   1e-5',...
            ' -scale_factor       0.04',...
            ' -adaptive           1',...
            ' ' opt];
    case {'canny_edges','break_edges','lines'}
        opt = [' -sigma              1',...
            ' -max_width          50',...
            ' -gauss_tail         0.0001',...
            ' -low                2',...
            ' -high               12',...
            ' -edge_min           60',...
            ' -min_length         10',...
            ' -border_size        2',...
            ' -border_value       0.0',...
            ' -scale              5.0',...
            ' -follow_strategy    2',...
            ' -join_flag          1',...
            ' -junction_option    0',...
            ' -bk_thresh          0.3',...
            ' -str_high           12',...
            ' -str_edge_min       60',...
            ' -str_min_length     60',...
            ' -min_fit_length     10',...
            ' ' opt];
    otherwise
        error('Unknown action required');
end

% Call the binary executable
fpath = fileparts(which(mfilename));
if strncmp(computer,'PC',2) % MS Windows
    exec_str = ['"' fpath '/xcv_segment.exe"'];
elseif strcmp(computer,'GLNX86') % Linux
    exec_str = [fpath '/xcv_segment_32bit'];
elseif strcmp(computer,'GLNXA64') % Linux
    exec_str = [fpath '/xcv_segment_64bit'];
else error('This function can run only with MS Windows or Linux');
end
%orig_dir = cd;
%cd(fileparts(which(mfilename)));
result = unix([exec_str  ' -i ' inname ' -f ' outname ' -' action ' ' opt]);
%cd(orig_dir);
if result ~= 0
    error('Calling the binary failed.');
end
try delete(inname); end

% Load the output file
f = fopen(outname,'r');
if f==-1
    error('Cannot load results from binary executable');
end
x = fscanf(f,'%g');
fclose(f);
try delete(outname); end

% Parse the output: see documentation to xcv_segment
c = 1;
switch action
    case 'corners'
        s = reshape( x(c+[1:2*x(1)]), [2 x(1)] );
        s = s + 1; % shift to matlab's (1,1)-based frame
        c = c + 1 + 2*x(1);
    case {'canny_edges','break_edges'}
        c = c + 1;
        for n = 1:x(1)
            s{n} = reshape( x(c+[1:4*x(c)]), [4 x(c)] );
            s{n}(1:2,:) = s{n}(1:2,:) + 1; % shift to matlab's (1,1)-based frame
            c = c + 1 + 4*x(c);
        end
    case 'lines'
        s = reshape( x(c+[1:4*x(1)]), [4 x(1)] );
        s = s + 1; % shift to matlab's (1,1)-based frame
        c = c + 1 + 4*x(1);
end
if length(x)+1 ~= c
    error('Wrong length of the output file');
end

return
