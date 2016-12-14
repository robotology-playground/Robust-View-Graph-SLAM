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
function [] = GenAveFeaSup(batchNumber, Nei, AbsFeaType, AbsFeaDate)

% This function calculate the feature of each subsuperpixel using texture
% infomation

if nargin < 1
	batchNumber = 1;
end

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
     ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
     Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName;

% ==================================================================================
% may chage for different usage
%load([ScratchDataFolder '/data/LowResImgIndexSuperpixelSep.mat']); % superpixel_index
%load([ScratchDataFolder '/data/DiffLowResImgIndexSuperpixelSep.mat']);
%load([ScratchDataFolder '/data/TextLowResImgIndexSuperpixelSep.mat']); 
% ==================================================================================
% load FeaMax to do normalizeing
load([GeneralDataFolder '/FeaMax.mat']);
FeaMax = 10.^floor(log10(FeaMax));

%load([ScratchDataFolder '/data/MaskGSky.mat']); % load maskg maskSky from CMU's output
%clear maskg;

% prepare data step
nu_pics = size(filename,2); % number of picture

batchImg = 1:batchSize:nu_pics;
l = 1;
for i = batchImg(batchNumber):min(batchImg(batchNumber)+batchSize-1, nu_pics)
%for i = [10:18 54:62]
       tic
       i 
% ==================================================================================
% may chage for different usage
%load([ScratchDataFolder '/data/LowResImgIndexSuperpixelSep.mat']); % superpixel_index
load([ScratchDataFolder '/data/DiffLowResImgIndexSuperpixelSep.mat']);
%load([ScratchDataFolder '/data/TextLowResImgIndexSuperpixelSep.mat']); 
% ==================================================================================
%       load([ScratchDataFolder '/data/MedSeg/MediResImgIndexSuperpixelSep' num2str(i) '.mat']);
%       Sup  = LowResImgIndexSuperpixelSep{i}; clear LowResImgIndexSuperpixelSep;
%       MedSup = MediResImgIndexSuperpixelSep; clear MediResImgIndexSuperpixelSep;
       DiffSup = DiffLowResImgIndexSuperpixelSep(i,end); clear DiffLowResImgIndexSuperpixelSep;
%       TextSup = TextLowResImgIndexSuperpixelSep(i,:,2:end); clear TextLowResImgIndexSuperpixelSep;
batchNumber,%       MedSup = imresize(MedSup, [vertical_size_hi_res horizontal_size_hi_res]);

       % load MediResImgIndexSuperpixelSep
       % decide to claen the Sup or not (means imclosing)
       %[Sup,MedSup]=CleanedSup(Sup,MedSup,maskSky{i});
       load([ScratchDataFolder '/data/CleanSup/CleanSup' num2str(i) '.mat']);
       % check if the MedSup and Sup have the same index
       % !!!!!!!Don't need to check it since later working on Sup scale only
%       Res = setdiff(unique(MedSup(:)),unique(Sup(:)));
%       if ~isempty(Res)
%           disp('error index from MedSup to Sup');
%           return;
%       end
       
       MedSup = double(MedSup);
       Sup = double(Sup);
       maskSky{i} = Sup == 0;  % the new skymap

       % load picsinfo just for the horizontal value
       PicsinfoName = strrep(filename{l},'img','picsinfo');
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

       % 1) Analyze Sup : SupFact with NuPatchEachSup and all 13 features calcuate in gen_fsup_new
       [SupFact, nList] = AnalyzeSup(Sup,maskSky{i}); % Calculate for H2
%       load([ScratchDataFolder '/data/temp/List' num2str(i) '.mat']);
       % ======================
      BounaryPHori = conv2(Sup,[1 -1],'same') ~=0;
      BounaryPHori(:,end) = 0;
      BounaryPVert = conv2(Sup,[1; -1],'same') ~=0;
      BounaryPVert(end,:) = 0;
      ClosestNList = [ Sup(find(BounaryPHori==1)) Sup(find(BounaryPHori==1)+VertYNuDepth);...
                       Sup(find(BounaryPVert==1)) Sup(find(BounaryPVert==1)+1)];
      ClosestNList = sort(ClosestNList,2);
      ClosestNList = unique(ClosestNList,'rows');
      ClosestNList(ClosestNList(:,1) == 0,:) = [];
%      MedBounaryPHori = conv2(MedSup,[1 -1],'same') ~=0;
%      MedBounaryPHori(:,end) = 0;
%      MedBounaryPVert = conv2(MedSup,[1; -1],'same') ~=0;
%      MedBounaryPVert(end,:) = 0;
%      MedClosestNList = [ MedSup(find(MedBounaryPHori==1)) MedSup(find(MedBounaryPHori==1)+SegVertYSize);...
%                          MedSup(find(MedBounaryPVert==1)) MedSup(find(MedBounaryPVert==1)+1)];
%      MedClosestNList = sort(MedClosestNList,2);
%      MedClosestNList = unique(MedClosestNList,'rows');
%      MedClosestNList(MedClosestNList(:,1) == 0,:) = [];
%       nList = [ClosestNList; MedClosestNList];
%       nList = unique(nList,'rows'); 
       nList = ClosestNList;
% =======================================================================================================
       % calculate all the features, ray, plane parameter, and row column value
       img = imread([GeneralDataFolder '/' ImgFolder '/' filename{l} '.jpg']);% read hi resolution image
       
       % check if the images resolusion is smaller then a certain size
       if prod(size(img))<SegVertYSize*SegHoriXSize*3
           disp('imresize hard work');
           img = imresize(img,[SegVertYSize SegHoriXSize],'bilinear');
       end    
       [vertical_size_hi_res horizontal_size_hi_res t] = size(img); clear t; 

       % generate lineseg
%       seglist=edgeSegDetection(img,i);       

       % generate the texture features 
       disp('going to cal Fea')

       % Start Sup relation with the 7 * 2 = 14 Multi Sup
       [MultiScaleSupTable] = MultiScalAnalyze( Sup, permute(  cat( 3, DiffSup{1,1}),...% DiffSup{1,2},...
                              [3 1 2])); 
%                              TextSup{1,1,1}, TextSup{1,1,2},...
%                              TextSup{1,2,1}, TextSup{1,2,2},...
%                              TextSup{1,3,1}, TextSup{1,3,2},...
%                              TextSup{1,4,1}, TextSup{1,4,2},...
%                              TextSup{1,5,1}, TextSup{1,5,2},...
%                              TextSup{1,6,1}, TextSup{1,6,2}),...
       clear DiffSup;% TextSup;

       global H2;
       [H2] = calculateFilterBanks_old(img); % (hard work 1min) use Ashutaosh's code
       clear img;
       H2 = permute(H2,[3 1 2]);

       % run Rajiv Code
%       FeaNList = findBoundaryFeaturesMore( Sup, MedSup, nList(:,1:2), seglist, FeaMax, 0);       
%       clear seglist;
 
       % Calculate the AveSupFea % H1 first
       [SupFea ] = AveSupFea(Sup, MedSup, SupFact, MultiScaleSupTable, FeaMax, 1); % Maight can be faster when doing calcuating fea from MultiScaleSupTable 
       H2 = H2.^2; % then H2
       [TempSupFea ] = AveSupFea(Sup, MedSup, SupFact, MultiScaleSupTable, FeaMax, 2); % Maight can be faster when doing calcuating fea from MultiScaleSupTable 
       SupFea = [SupFea TempSupFea(:,2:end)]; % Maight can be faster when doing calcuating fea from MultiScaleSupTable 
       H2 = H2.^2; % then H4
       TempSupFea = AveSupFea(Sup, MedSup, SupFact, MultiScaleSupTable, FeaMax, 4);
       SupFea = [SupFea TempSupFea(:,2:end)]; % Maight can be faster when doing calcuating fea from MultiScaleSupTable 
       clear MultiScaleSupTable TempSupFea;
       clear global H2;
       toc

       % gnerate Feature as the same row of nList , and 
       tic
       if Nei
          load([ScratchDataFolder '/data/feature_Abs_' AbsFeaType int2str(batchNumber) '_' AbsFeaDate '.mat']); % 'f'
       else
          f = [];
       end

       % add new option to include Col feature 103*2 col features and 103*2 row features and 103*4 Nei features 
       [FeaNList, nList] = GenFeaParaNList( Sup, MedSup, maskSky{i}, nList, SupFact, SupFea, FeaMax, Nei, f{i-10*(batchNumber-1)}, i);
       clear f;
%       FeaNList = [ FeaNListTemp FeaNList];
       toc  


       disp([ScratchDataFolder '/data/SupFea/FeaNList' num2str(i) 'new.mat']); 
       save([ScratchDataFolder '/data/SupFea/FeaNList' num2str(i) 'new.mat'],'FeaNList','nList');
%return;
%       save([ScratchDataFolder '/data/SupFea/DiffA_Alpha' num2str(i) '.mat'],'DiffA','DiffAlpha');
       clear FeaNListTemp nList Sup MedSup nList SupFact SupFea FeaNList;  
       l= l+1;       
end
%save([ScratchDataFolder '/data/SupFea/FeaNList' num2str(i) '.mat'],'FeaNList','nList');
%save([ScratchDataFolder '/data/SupFea/DiffA_Alpha' num2str(i) '.mat'],'DiffA','DiffAlpha');
return;
