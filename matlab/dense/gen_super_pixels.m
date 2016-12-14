
function [MedSup, Sup, param, SupNeighborTableFlag] = gen_super_pixels(param, img, SelectSegmentationPara)
% this function generate superpixel using default parameter
% but can also change to manually input parameter

%%% Jeff's Comments (Min modified)
% Send that image to the CMU
% segmentation program with params 0.8*[sigm, k, min].  If
% SelectSegmenationPara is true, then display the results and ask
% the user to enter new sigm, k, and min; repeat until user is happy.
%
% Output the MedSup(only the smallest scale) and Sup(three scale)
%%%%

% % Input parameters
% param.display : show superpixels ?
% param.SegVertYSize (900) :
% param.SegHoriXSize (1200) :
% param.sigm (.5, .8, .3) :
% param.k (100. 200, 300) :
% param.minp (100, 150, 20) :
% param.OutPutFolder :
% param.filename :
% param.PpmOption (1) : option to storage the ppm image segmentation
% param.VertYNuDepth (55) :
% param.HoriXNuDepth (305) :
% img :
% SelectSegmentationPara : if SelectSegmentationPara == 1, enable the parameter interation with user.
            
% segment code path
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/dense/superpixels/code/segment');

% default parameter
if nargin < 3
    SelectSegmentationPara = 0;
end
display = param.display; % set to display or not

% use different scale to generate small(0.8) middle(1.6) 5(large) scale of superpixel
scale = param.scale;

% find the dimension size of the Hi Resolution image
[VertYSizeHiREs, HoriXSizeHiREs, ~]= size(img);

% using a fixed range of median size image [SegVertYSize SegHoriXSize]
%  to generate superpixel to reduce computation
if VertYSizeHiREs*HoriXSizeHiREs > param.SegVertYSize*param.SegHoriXSize
    % downsample high resolution image to a fixed median size image
    % +4 because segmentImgOpt gives constant additinal rows and column
    % so add 4 rows and columns a prior then subtract them
    img = imresize(img,[param.SegVertYSize+4 param.SegHoriXSize+4 ],'nearest');
else
    param.SegVertYSize = VertYSizeHiREs-4;
    param.SegHoriXSize = HoriXSizeHiREs-4;
end

% generate superpixel of each image
for j = 1:3% number of scale of superpixel
    
    % choose superpixel of the images
    % default segmentation parameter
    ok = 0; % ok ==1 means accept the segmentation
    while 1
        
        % call the efficient segment function writen in C++ from MIT
        % Output the high resolution image ( + 1 since the smallest index can be zero)
        if j ==1
            a = segmentImgOpt( param.sigm*scale(j), param.k*scale(j), param.minp*scale(j), img,...
                [ param.OutPutFolder param.filename '.ppm'], param.PpmOption) + 1;
        else
            a = segmentImgOpt( param.sigm*scale(j), param.k*scale(j), param.minp*scale(j), img,...
                [ param.OutPutFolder param.filename '.ppm'], 0) + 1;
        end
        a = a(3:(end-2),3:(end-2)); % clean the edge superpixel index errors
        
        % arrange the superpixel index in order
        if j == 1 % For the smallest Scale
            
            ma = max(a(:));
            Unique_a = unique(a);
            SparseIndex = sparse(ma,1);
            SparseIndex(Unique_a) = 1:size(Unique_a);
            MedSup = full(SparseIndex(a));
            
            % downsample to size as prediected depth map
            Sup{j} = imresize(MedSup,[param.VertYNuDepth param.HoriXNuDepth],'nearest');
            % clean superpixel section
            % merage all small and disconneted points in 1st scale segmentation
            [Sup{j}, SupNeighborTableFlag] = merge_superpixels(Sup{j}, param);
            
        else  % o/w don't need the MedSup
            
            % downsample to size size as prediected depth map
            a = imresize(a,[param.VertYNuDepth param.HoriXNuDepth],'nearest');
            ma = max(a(:));
            Unique_a = unique(a);
            SparseIndex = sparse(ma,1);
            SparseIndex(Unique_a) = 1:size(Unique_a);
            Sup{j} = full(SparseIndex(a));
        end
        clear a SparseIndex Unique_a ma;
        
        
        % show superpixel
        if display
            figure(1);
            imagesc(Sup{j});
            newmap = rand(max(max(Sup{j})),3);
            colormap(newmap);
        end
        
        % check if need to select segmentation parameter
        if SelectSegmentationPara==1;
            ok = input('Is the segmentation of image OK');% input new segmentation parameter
        else
            ok =1 ;% accept default segmentation parameter
        end
        
        if ok==1;
            break;
        end
        
        % Get the user selected parameter
        sigm = input('type sigm of segmentation');
        k = input('type k of segmentation');
        minp = input('type min of segmentation');
        
    end % end of while 1
    
end % end of for j=1:3

return;
