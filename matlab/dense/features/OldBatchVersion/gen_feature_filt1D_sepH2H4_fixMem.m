% *  This code was used in the following articles:
% *  [1] Learning 3-D Scene Structure from a Single Still Image, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on 3D Representation for Recognition (3dRR-07), 2007.
% *      (best paper)
% *  [2] 3-D Reconstruction from Sparse Views using Monocular Vision, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on Virtual Representations and Modeling 
% *      of Large-scale environments (VRML), 2007. 
% *  [3] 3-D Depth Reconstruction from a Single Still Image, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      International Journal of Computer Vision (IJCV), Aug 2007. 
% *  [6] Learning Depth from Single Monocular Images, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      In Neural Information Processing Systems (NIPS) 18, 2005.
% *
% *  These articles are available at:
% *  http://make3d.stanford.edu/publications
% * 
% *  We request that you cite the papers [1], [3] and [6] in any of
% *  your reports that uses this code. 
% *  Further, if you use the code in image3dstiching/ (multiple image version),
% *  then please cite [2].
% *  
% *  If you use the code in third_party/, then PLEASE CITE and follow the
% *  LICENSE OF THE CORRESPONDING THIRD PARTY CODE.
% *
% *  Finally, this code is for non-commercial use only.  For further 
% *  information and to obtain a copy of the license, see 
% *
% *  http://make3d.stanford.edu/publications/code
% *
% *  Also, the software distributed under the License is distributed on an 
% * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
% *  express or implied.   See the License for the specific language governing 
% *  permissions and limitations under the License.
% *
% */
function []=gen_feature_filt1D_sepH2H4_fixMem(batchNumber,HistFeaType,Absolute)

% This function calculate the feature of each subsuperpixel using texture
% infomation

% decide the Hist
if strcmp(HistFeaType,'Whole')
   Hist = 1;
   HistFeaType
else
   Hist = 0;
end

if Hist ~= 1 && Absolute ~=1
	return
end

if nargin < 1
	batchNumber = 1;
elseif nargin < 2
	Hist = 1;
        Absolute =1;
elseif nargin < 3 
        Absolute =1;
end

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName...
    TrainVerYSize TrainHoriXSize MempryFactor;


% load estimated sky
load([ScratchDataFolder '/data/MaskGSky.mat']); % maskg is the estimated ground maskSky is the estimated sky

% load([ScratchDataFolder '/data/filename.mat']);% load the filename 
% load([ScratchDataFolder '/data/PlaneParameterTure.mat']); % planeParameter 
load([ScratchDataFolder '/data/LowResImgIndexSuperpixelSep.mat']); % superpixel_index
%load([ScratchDataFolder '/data/MediResImgIndexSuperpixelSep.mat']); % MediResImgIndexSuperpixelSep
load([ScratchDataFolder '/data/DiffLowResImgIndexSuperpixelSep.mat']); % DiffLowResImgIndexSuperpixelSep
DiffLowResImgIndexSuperpixelSep = DiffLowResImgIndexSuperpixelSep(:,1);% need only the middle scale segmentation

%load([ScratchDataFolder '/data/FeatureSuperpixel.mat']); %load feature of superpixel
% load FeaMax to do normalizeing
load([GeneralDataFolder '/FeaMax.mat']);
FeaMax = 10.^floor(log10(FeaMax));

% prepare data step
nu_pics = size(filename,2); % number of pictures


% ====================== change able parameter =====================
%batchSize = 10;% decide the batch size as 10 image pre batch
% ==================================================================

% for batchImg = 1:batchSize:nu_pics
batchImg = 1:batchSize:nu_pics;

f = []; % total feature: feature of superpixel followed by texture feature of patch
    PicsNu = 1
    for i = batchImg(batchNumber):min(batchImg(batchNumber)+batchSize-1, nu_pics)
       i 
       % load MediResImgIndexSuperpixelSep
       load([ScratchDataFolder '/data/MedSeg/MediResImgIndexSuperpixelSep' num2str(i) '.mat']);        
  
       % calculate all the features, ray, plane parameter, and row column value
       img = imread([GeneralDataFolder '/' ImgFolder '/' filename{i} '.jpg']);% read in the hi resolution image of the ith file
       size(img) 
       % change the resolution in to exactly [TrainVerYSize TrainHoriXSize]
       
       if ~all(size(img) == [TrainVerYSize TrainHoriXSize 3])
           disp('resize to 2272 1704')
           img = imresize(img,[TrainVerYSize TrainHoriXSize],'bilinear');
       end
       % the images resolusion must be bigger then a certain size to have reasonable predicted depth
%       if prod(size(img))<SegVertYSize*SegHoriXSize*3
%           img = imresize(img,[SegVertYSize SegHoriXSize],'bilinear');
       % the images resolusion must be smaller then a certain size to avoid out of memory
%       elseif prod(size(img)) > TrainVerYSize*TrainHoriXSize*3*(MempryFactor);
%           disp('origin size')
%           size(img)
%           img = imresize(img,[TrainVerYSize TrainHoriXSize],'bilinear');
%           disp('image too big')
%           size(img)
%       end    
       [vertical_size_hi_res horizontal_size_hi_res t] = size(img); clear t; 
	   % get the horizontal(vertical_size_hi_res) and vertical(horizontal_size_hi_res) size of hi resolution image
       sup_hi_res = imresize(MediResImgIndexSuperpixelSep, ...
       [vertical_size_hi_res horizontal_size_hi_res],'nearest');% enlarge the low res superpixel into hi res superpixel
       clear MediResImgIndexSuperpixelSep;  
 
       % generate the superpixel in the depth_grid size
       sup = imresize(LowResImgIndexSuperpixelSep{i,1},[VertYNuDepth HoriXNuDepth],'nearest');       
       
       % load picsinfo just for the horizontal value
        PicsinfoName = strrep(filename{i},'img','picsinfo');
        temp = dir([GeneralDataFolder '/PicsInfo/' PicsinfoName '.mat']);
        if size(temp,1) == 0
            a = a_default;
            b = b_default;
            Ox = Ox_default;
            Oy = Oy_default;
            Horizon = Horizon_default;
        else    
            load([GeneralDataFolder '/PicsInfo/' PicsinfoName '.mat']);
        end
        
       % calculate the ray 
       RayCenter = GenerateRay(HoriXNuDepth,VertYNuDepth,'center',a,b,Ox,Oy); %[ horiXSizeLowREs VertYSizeLowREs 3]
        
       % calculate how many mask we need
       NuMask = ceil([HoriXNuDepth/HoriXNuPatch; VertYNuDepth/VertYNuPatch])

       % Grid Info
       gridinfo = [horizontal_size_hi_res HoriXNuDepth HoriXNuPatch;	vertical_size_hi_res VertYNuDepth VertYNuPatch];
       ratio(1:2) = gridinfo(:,1)./gridinfo(:,end);
       ratio(3:4) = gridinfo(:,1)./gridinfo(:,2);

       % calcuate the position of the mask
       hight(1) = round((floor(ratio(2))-1)/2);
       hight(2) = floor(ratio(2)) - 1 - hight(1);
       width(1) = round((floor(ratio(1))-1)/2);
       width(2) = floor(ratio(1)) - 1 - width(1);

       %big_sup_depthmap_res = imresize(LowResImgIndexSuperpixelSep{i,2}, ...
	%		[VertYNuDepth HoriXNuDepth],'nearest');% enlarge the low res superpixel into hi res superpixel
       f_pics = [];
       tic
       [H2] = calculateFilterBanks_old(img); % (hard work 1min) use Ashutaosh's code
       H2 = H2.^2;
%        NuFeaH2H4 = 2*size(H2,3);
       hcol = ones(floor(ratio(2)),1);
       hrow = ones(1,floor(ratio(1)));
       hrect = ones(floor(ratio(2)),floor(ratio(1)));
size(hcol)
size(hrow)
size(hrect)

if Absolute == 1
       % highest resolution
       row_start = 1;
%        MinTest = zeros(55*305,1);
       for j = 1:NuMask(1)
           for k = 1:NuMask(2);
               f_pics_mask = [];
               % first generate the mask respect to the dominate
               % subsuperpixel 
               [mask,SupIndex,PixelMask,PatchMask] = makeSubSupMask(gridinfo,sup_hi_res,sup,[j; k],width,hight);
               % calculaing the normalize value
                  NormalizeValue = conv2(hcol,hrow,mask,'same');

               % 1) zeroth feature the superpixel index to keep a record
               f_pics_mask = [f_pics_mask SupIndex];
               clear SupIndex;

               % 2) generate the 1:34 features for H2  for 1 center and 4 neighbor (left right top bottom)
               fInd=2;
               for m = 1:17
                   %temp1 = conv2(hcol,hrow,H2(:,:,m).*(mask),'same');
                   temp = conv2(H2(:,:,m).*(mask),hrect,'same');
                   size(temp)
                   f_pics_mask = [f_pics_mask temp(PixelMask)./NormalizeValue(PixelMask)...
                                  ./FeaMax(1,fInd)];
                   fInd = fInd+1;
               end   

               f_pics(PatchMask, row_start:row_start+size(f_pics_mask,2)-1) = f_pics_mask; 
               
           end
       end

       H2 = H2.^2; % This H2 is H4
       row_start = row_start+size(f_pics_mask,2);

	% RAJIV MIN -- ERROR ------- this part which sums up the filter response in a superpixel, should
	% be done using integral images, prefably in C++ code.

       for j = 1:NuMask(1)
           for k = 1:NuMask(2);
               f_pics_mask = [];
               % first generate the mask respect to the dominate
               % subsuperpixel 
               [mask,SupIndex,PixelMask,PatchMask] = makeSubSupMask(gridinfo,sup_hi_res,sup,[j; k],width,hight);
               % calculaing the normalize value
                  NormalizeValue = conv2(hcol,hrow,mask,'same');

               % 1) zeroth feature the superpixel index to keep a record
%               f_pics_mask = [f_pics_mask SupIndex];
               clear SupIndex;

               % 2) generate the 1:34 features for H4 for 1 center and 4 neighbor (left right top bottom)
               fInd=19;
               for m = 1:17
                   %temp = conv2(hcol,hrow,H2(:,:,m).*(mask),'same');
                   temp = conv2(H2(:,:,m).*(mask),hrect,'same');
                   f_pics_mask = [f_pics_mask temp(PixelMask)./NormalizeValue(PixelMask)...
                                 ./FeaMax(1,fInd)];
                   fInd = fInd+1;
               end

               f_pics(PatchMask, row_start:row_start+size(f_pics_mask,2)-1) = f_pics_mask; 
               
           end
       end
end

if Hist == 1
% RAJIV ---- ERROR -- should be disabled.
H2 = H2.^(0.5);
% =============calculate the histagram of the features for relative depth estimation================
        disp('cal_relative')
	[relativeFeatureVector] = makeRelativeFeatureVector(H2,1);
% ==================================================================================================
end
       clear H2;

       % 1/3 resolution
       feaScale1 = [];
       MedResY = round(gridinfo(2,1)/3);      
       MedResX = round(gridinfo(1,1)/3);
       DepthGridSizeY = MedResY/VertYNuDepth;      
       DepthGridSizeX= MedResX/HoriXNuDepth;      
       imgMedRes = imresize(img,[MedResY MedResX],'nearest');
       clear img;
       [H2] = calculateFilterBanks_old(imgMedRes); % (hard work 1min) use Ashutaosh's code
       H2 = H2.^2;
       H4 = H2.^2;
if Absolute == 1
       row_start = size(f_pics,2)+1;
       % calculating number of mask 
       NormalizeValue = conv2(hcol,hrow,ones(MedResY,MedResX),'same');
       NuMask = ceil([HoriXNuDepth/HoriXNuPatch*3; VertYNuDepth/VertYNuPatch*3]); 
              
       % 1) generating the PixelMask
       PixelMask = logical(zeros(MedResY,MedResX));
       [X Y] = meshgrid(ceil((1/2)*DepthGridSizeX:DepthGridSizeX:MedResX),...
                        ceil((1/2)*DepthGridSizeY:DepthGridSizeY:MedResY));
       PixelMask = sub2ind(size(PixelMask),Y(:),X(:));
       for m = 1:17
           temp = conv2(hcol,hrow,H2(:,:,m),'same');% 102.533776 seconds
           feaScale1(:,m) = temp(PixelMask)./NormalizeValue(PixelMask)...
               ./FeaMax(1,fInd);
           fInd = fInd +1;
       end   
               
       for m = 1:17
           temp = conv2(hcol,hrow,H4(:,:,m),'same');% 102.533776 seconds
           feaScale1(:,m+17) = temp(PixelMask)./NormalizeValue(PixelMask)...
               ./FeaMax(1,fInd);
           fInd = fInd +1;
       end

%       shift = [0 0; -1 0; 1 0; 0 -1; 0 1].*repmat(NuMask',[5 1]);
%       for l = 1:5
%           [Ix Iy] = meshgrid(max(min(2+shift(l,1):HoriXNuDepth+2-1+shift(l,1),HoriXNuDepth+2),1),...
%                     max(min(2+shift(l,2):VertYNuDepth+2-1+shift(l,2),VertYNuDepth+2),1));
%           maskNeibor = sub2ind([VertYNuDepth+2, HoriXNuDepth+2], Iy(:), Ix(:));
%           f_pics = [f_pics  feaScale1(maskNeibor,:)];
%       end
       
       size(f_pics)
       size(feaScale1)    
       f_pics = [f_pics  feaScale1];
end
clear H4;

if Hist == 1
% =============calculate the histagram of the features for relative depth estimation================
	[relativeFeatureVector] = cat(3,relativeFeatureVector,makeRelativeFeatureVector(H2,2));
% ==================================================================================================
end

       clear H2;

       % 1/9 resolution
       feaScale1 = [];
       LowResY = round(gridinfo(2,1)/9);      
       LowResX = round(gridinfo(1,1)/9);      
       DepthGridSizeY = LowResY/VertYNuDepth;      
       DepthGridSizeX= LowResX/HoriXNuDepth;      
       imgLowRes = imresize(imgMedRes,[LowResY LowResX],'nearest');
       clear imgMedRes;
       [H2] = calculateFilterBanks_old(imgLowRes); % (hard work 1min) use Ashutaosh's code
       H2 = H2.^2;
       H4 = H2.^2;

if Absolute == 1
       row_start = size(f_pics,2)+1;
       % calculating number of mask 
       NormalizeValue = conv2(hcol,hrow,ones(LowResY,LowResX),'same');
       NuMask = ceil([HoriXNuDepth/HoriXNuPatch*9; VertYNuDepth/VertYNuPatch*9]); 
               
       % 1) generating the PixelMask
       PixelMask = logical(zeros(LowResY,LowResX));
       [X Y] = meshgrid(ceil((1/2)*DepthGridSizeX:DepthGridSizeX:LowResX),...
                        ceil((1/2)*DepthGridSizeY:DepthGridSizeY:LowResY));
       PixelMask = sub2ind(size(PixelMask),Y(:),X(:));
%       size(maskInd)
%       PixelMask(maskInd) = true;
       %PixelMask = logical(zeros(LowResY,LowResX));
       %PixelMask(ceil((1/2)*DepthGridSizeY:DepthGridSizeY:LowResY),...
       %          ceil((1/2)*DepthGridSizeX:DepthGridSizeX:LowResX)) = true;
%       valuemask = sum(PixelMask,1)>0;
%       if all(PixelMask(1,valuemask) == true)
%          PixelMask(2,valuemask) = true;
%       else
%          PixelMask(1,valuemask) = true;
%       end   
%       if all(PixelMask(end,valuemask) == true)
%          PixelMask(end-1,valuemask) = true;
%       else
%          PixelMask(end,valuemask) = true;
%       end   
%          valuemask = sum(PixelMask,2)>0;
%       if all(PixelMask(valuemask,1) == true)
%          PixelMask(valuemask,2) = true;
%       else
%          PixelMask(valuemask,1) = true;
%       end   
%       if all(PixelMask(valuemask,end) == true);
%          PixelMask(valuemask,end-1) = true;
%       else
%          PixelMask(valuemask,end) = true;
%       end               
               
       % 2) generate the 1:34 features for H2 and H4 for 1 center and 4 neighbor (left right top bottom)
       for m = 1:17
           temp = conv2(hcol,hrow,H2(:,:,m),'same');% 102.533776 seconds
           feaScale1(:,m) = temp(PixelMask)./NormalizeValue(PixelMask)...
                   ./FeaMax(1,fInd);
           fInd = fInd +1;
       end   
               
       for m = 1:17
           temp = conv2(hcol,hrow,H4(:,:,m),'same');% 102.533776 seconds
           feaScale1(:,m+17) = temp(PixelMask)./NormalizeValue(PixelMask)...
                   ./FeaMax(1,fInd);
           fInd = fInd +1;
       end

%       shift = [0 0; -1 0; 1 0; 0 -1; 0 1].*repmat(NuMask',[5 1]);
%       for l = 1:5
%           [Ix Iy] = meshgrid(max(min(2+shift(l,1):HoriXNuDepth+2-1+shift(l,1),HoriXNuDepth+2),1),...
%%                     max(min(2+shift(l,2):VertYNuDepth+2-1+shift(l,2),VertYNuDepth+2),1));
%           maskNeibor = sub2ind([VertYNuDepth+2, HoriXNuDepth+2], Iy(:), Ix(:));
%           f_pics = [f_pics  feaScale1(maskNeibor,:)];
%       end
       f_pics = [f_pics  feaScale1];
end          
clear H4;

if Hist ==1
% =============calculate the histagram of the features for relative depth estimation================
	[relativeFeatureVector] = cat(3,relativeFeatureVector,makeRelativeFeatureVector(H2,5));
        RFVector{PicsNu} = relativeFeatureVector;
        clear relativeFeatureVector;
        DateStamp = date;
        save([ScratchDataFolder '/data/feature_Hist_Whole' num2str(batchNumber) '_' DateStamp  '.mat'],'RFVector');
% ==================================================================================================
end

       clear H2;
 
if Absolute == 1
       % superpixel features
       %fsup = FeatureSuperpixel{i};
       %f_pics = [f_pics fsup(:,f_pics(:,1))'];
       % other features without relation with neiborfeatures
       % 1) the closest ground position to the (i, j) patch
         % How can we remove DiffLowResImgIndexSuperpixelSep since generate another superpixel takes time
         big_sup = imresize(DiffLowResImgIndexSuperpixelSep{i,1},[VertYNuDepth HoriXNuDepth],'nearest');
         GroundSupIndex = unique(big_sup(maskg{i}));

%         NonGround = [];
%       for j = 1:HoriXNuDepth
%           Ground(j) = analysesupinpatch(big_sup(round(gridinfo(2,2)*(1-Horizon)):end,j));
%            NonGround(j) = analysesupinpatch(big_sup(round(1:(gridinfo(2,2)*(1-Horizon)-1)),j));
%           NonGround = [NonGround (unique(big_sup(round(1:(gridinfo(2,2)*(1-Horizon)-1)),j)))'];
%       end
%       Ground = setdiff(Ground,NonGround);
       for j = 1:HoriXNuDepth
           Gmask = logical(zeros(size( big_sup(:,j))));
           for k=GroundSupIndex'
               Gmask(big_sup(:,j) == k) = true;
           end
           [rowSub colSub ] = find(Gmask);
           minRow = min(rowSub);
%           if size(minRow,1) == 0
%               Gmask = logical(zeros(size( big_sup(:,j))));
%               MaybeGround = unique(big_sup(round(gridinfo(2,2)*(1-Horizon)):end,j));
%               MaybeGround = setdiff(MaybeGround,NonGround);
%               for k=MaybeGround
%                   Gmask(big_sup(:,j) == k) = true;
%               end
%               [rowSub colSub ] = find(Gmask);
%               minRow = min(rowSub);
%               if size(minRow,1) == 0
%                  lastguessGround = analysesupinpatch(big_sup(round(gridinfo(2,2)*(1-Horizon)):end,j));
%                  [rowSub colSub ] = find(big_sup(round(gridinfo(2,2)*(1-Horizon)):end,j)==lastguessGround);
%                  minRow = min(rowSub);
%               end
%           end    
           %minRow = max(minRow,round(gridinfo(2,2)*(1-Horizon)));

           GroundVertEdge(:,j) = (1-(VertYNuDepth-1)/VertYNuDepth)./RayCenter(:,j,3);
           if size(minRow,1) ~= 0
           	GroundVertEdge((minRow+1):VertYNuDepth,j) = (1-((minRow+1):VertYNuDepth)'/VertYNuDepth)./RayCenter((minRow+1):VertYNuDepth,j,3);
                GroundVertEdge(1:(minRow),j) = (1- minRow/VertYNuDepth)./RayCenter(1:(minRow),j,3);
           end
       end
       GroundVertEdgePics = GroundVertEdge;


       %f_pics(:,end+1) = max([GroundVertEdgePics(:) repmat((1:VertYNuDepth)',[HoriXNuDepth 1])],[],2);
       %f_pics(:,end+1) = Gr undVertEdgePics(:);%./reshape(RayCenter(:,:,3),[],1);
       f_pics(:,end+1) = min([GroundVertEdgePics(:)...
                repmat((VertYNuDepth:-1:1)'/VertYNuDepth,[HoriXNuDepth 1])./reshape(RayCenter(:,:,3),[],1)],[],2);    
       
       disp(['Image Number ' num2str(i)]);
       f{PicsNu} = f_pics;
       % Absolute Features are index by Type BatchNu Data
       DateStamp = date;
       save([ScratchDataFolder '/data/feature_Abs_Whole' num2str(batchNumber) '_' DateStamp  '.mat'],'f');

toc;        
%return;
end
       PicsNu = PicsNu + 1
       batchNumber
    end
    % Absolute Features are index by Type BatchNu Data
    DateStamp = date;
    if Absolute == 1
    save([ScratchDataFolder '/data/feature_Abs_Whole' num2str(batchNumber) '_' DateStamp  '.mat'],'f');
    end
    % Hist Features are index by Type BatchNu Data
    if Hist == 1
       save([ScratchDataFolder '/data/feature_Hist_Whole' num2str(batchNumber) '_' DateStamp  '.mat'],'RFVector');
    end
return;
