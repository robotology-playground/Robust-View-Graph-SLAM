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
function []=gen_feature_general(tskName,imgFolder,trainSet,learnType,learnSkyEx,learnLog,learnNear,...
                   learnAlg,learnDate,absFeaType,absFeaDate,HistFeaType,histFeaDate, ...
                   generalDataFolder, scratchDataFolder, localFolder, clusterExecutionDir, batchNumber,Absolute)

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
%load([ScratchDataFolder '/data/MaskGSky.mat']); % maskg is the estimated ground maskSky is the estimated sky

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
f_pics = [];
    PicsNu = 1
    for i = batchImg(batchNumber):min(batchImg(batchNumber)+batchSize-1, nu_pics)
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
     
        RayCenter = GenerateRay(HoriXNuDepth,VertYNuDepth,'center',a,b,Ox,Oy); %[ horiXSizeLowREs VertYSizeLowREs 3]

        Default=SetupDefault(tskName,imgFolder,trainSet,learnType,learnSkyEx,learnLog,learnNear,learnAlg,learnDate,absFeaType, ...
            absFeaDate,HistFeaType,histFeaDate,generalDataFolder,scratchDataFolder,localFolder,clusterExecutionDir); 
        [TextureFeature]=GenTextureFeature(Default, img, sup, sup_hi_res, 1);
        
        f_pics=TextureFeature.Abs;

       %% write ground boundary here
     if Absolute == 1
       % superpixel features
       %fsup = FeatureSuperpixel{i};
       %f_pics = [f_pics fsup(:,f_pics(:,1))'];
       % other features without relation with neiborfeatures
       % 1) the closest ground position to the (i, j) patch
         % How can we remove DiffLowResImgIndexSuperpixelSep since generate another superpixel takes time
         big_sup = imresize(DiffLowResImgIndexSuperpixelSep{i,1},[VertYNuDepth HoriXNuDepth],'nearest');
%%%%%         GroundSupIndex = unique(big_sup(maskg{i}));
         DefaultGroundMask = [zeros(floor(VertYNuDepth/2), HoriXNuDepth); ones(VertYNuDepth-floor(VertYNuDepth/2), HoriXNuDepth)];
         GroundSupIndex = unique(big_sup(logical(DefaultGroundMask)));

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
       %size(f_pics)
       %size(sup(:))
       f{PicsNu} = [sup(:) f_pics];

     end  

       PicsNu = PicsNu + 1
       batchNumber

    end 
    save([ScratchDataFolder '/data/feature_Abs_Whole' num2str(batchNumber) '_.mat'],'f');
return;
