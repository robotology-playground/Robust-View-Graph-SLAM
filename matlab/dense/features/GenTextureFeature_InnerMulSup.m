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
function [TextureFeature, TextSup]=GenTextureFeature_InnerMulSup(Default, Img, SupMedScale, SmallSup, HiSup, FeatureType)%, maskg)

% This function calculate all features that using the
% output from calculateFilterBanks_old (Time and ram comsuming)
% These Feature including AbsFeature DiffFeature HistogramFeature
% FeatureTypr:            1          2           3 

% Interface:
% Input---
% 1) Default is the structure contain all the default parameter setting
% 2) Img is the RGB image matrix
% 3) SmallSup is Superpixel index matrix of size Default.VertYNuDepth Default.HoriXNuDepth
% 4) HiSup is HiResolution Superpixel index matrix of size Default.SegVertYSize Default.SegHoriXSize or smaller
% 5) FeatureType AbsFeature(1) DiffFeature(2) HistogramFeature(3)
% Output---
% 1) TextureFeature.Abs:
% 2) TextureFeature.Diff:
% 3) TextureFeature.Hist:

% Flag set up for testing purpose ========================================
% Must be stablized after testing XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Default.Flag.NormalizeFlag = 1; % Default.Flags.NormalizeFlag set to 1 if we want to normalize the Feature according to the Trainning Set
DominateSupFilter = 1; % DominateSupFilter set to 1 if we want to use only the dominate superpixel area for acerage feature calculation
% ========================================================================

% load FeaMax to do normalizeing
load([Default.ParaFolder 'FeaMax.mat']);
if Default.Flag.NormalizeFlag == 1
   FeaMax = 10.^floor(log10(FeaMax));
else
   FeaMax(:) = 1;
end

% Initial Setting of TextureFeature Structure:
% size(TextureFeature.Abs) = (No of Depth Point, No of Features: 17 texture filters, H2 and H4 kinds, 3 Scales)
TextureFeature.Abs = zeros(Default.VertYNuDepth*Default.HoriXNuDepth, 17*2*3);
TextureFeature.Diff = zeros(1,1);
TextureFeature.Hist = zeros(1,1);

% change the resolution in to exactly [TrainVerYSize TrainHoriXSize]

if ~all(size(Img) == [Default.TrainVerYSize Default.TrainHoriXSize 3])
   if Default.Flag.DisplayFlag
      disp('Resize to 2272 1704');
   end
   Img = imresize(Img,[Default.TrainVerYSize Default.TrainHoriXSize],'bilinear');
end

% Start the loop for 3 Scale
fInd = 2;
for i = 1:3
    
    if DominateSupFilter && i==1 % only if DominateSupFilter and the largest scale will use superpixel mask
       SupMaskFlag = 1;
    else 
       SupMaskFlag = 0;
    end
    
    % Generate the 17 texture filter output
    % set to global for InnerMultiple sup generation
    global H2;
    [H2] = calculateFilterBanks_old(Img); % (hard work 1min) use Ashutaosh's code
    
    % InnerMulSup =================================================================
    if i == 2        
       [TextSup]=gen_TextSup_efficient_GlobalH2( Default);
    end
    % =====================================================================
    
    % Scale Setting
    ResY = round(Default.TrainVerYSize/(3^(i))); % ResY ResX will be the new image size     
    ResX = round(Default.TrainHoriXSize/(3^(i)));
    Img = imresize(Img,[ResY ResX],'nearest');

    H2 = H2.^2; % Now H2 stand for H2 22222222222222222222222222222222222222222222222222222222222222222222222222222
 
    if any(FeatureType == 1) % calculate the AbsFeatures
       StartAbs = 34*(i-1)+1;
%       [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)) fInd]= AbsFeatureGenMex(Default, H2, SmallSup, HiSup, SupMaskFlag, FeaMax, fInd);
       [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)) fInd]= ...
            AbsFeatureGenMex_MemoryEfficient(Default, SmallSup, HiSup, ...
            SupMaskFlag, FeaMax, fInd);
    end

    if any(FeatureType == 3) % calculate the HistogramFeature
       [relativeFeatureVector] = makeRelativeFeatureVector(H2,1);
    end

    if any(FeatureType == 2) % calculate the DiffFeature

    end

    % 2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222

    H2 = H2.^2; % This H2 is H4 44444444444444444444444444444444444444444444444444444444444444444444444444444444444

    if any(FeatureType == 1) % calculate the AbsFeatures
       StartAbs = 34*(i-1)+1+17;
%       [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)) fInd]= AbsFeatureGenMex(Default, H2, SmallSup, HiSup, SupMaskFlag, FeaMax, fInd);
       [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)) fInd]= AbsFeatureGenMex_MemoryEfficient(Default, SmallSup, HiSup, SupMaskFlag, FeaMax, fInd);
    end
    % 444444444444444444444444444444444444444444444444444444444444444444444
    
    clear global H2;
end

% GoundBoundaryFea
         RayCenter = GenerateRay(Default.HoriXNuDepth,Default.VertYNuDepth,'center',...
                     Default.a_default,Default.b_default,Default.Ox_default,...
                     Default.Oy_default); %[ horiXSizeLowREs VertYSizeLowREs 3]
        big_sup = imresize(SupMedScale,[Default.VertYNuDepth Default.HoriXNuDepth],'nearest');
        DefaultGroundMask = [zeros(floor(Default.VertYNuDepth/2), Default.HoriXNuDepth); ones(Default.VertYNuDepth-floor(Default.VertYNuDepth/2), Default.HoriXNuDepth)];
        GroundSupIndex = unique(big_sup(logical(DefaultGroundMask)));
%          GroundSupIndex = unique(big_sup(logical(maskg)));
       for j = 1:Default.HoriXNuDepth
           Gmask = logical(zeros(size( big_sup(:,j))));
           for k=GroundSupIndex'
               Gmask(big_sup(:,j) == k) = true;
           end
           [rowSub colSub ] = find(Gmask);
           minRow = min(rowSub);

           GroundVertEdge(:,j) = (1-(Default.VertYNuDepth-1)/Default.VertYNuDepth)./RayCenter(:,j,3);
           if size(minRow,1) ~= 0
                GroundVertEdge((minRow+1):Default.VertYNuDepth,j) = (1-((minRow+1):Default.VertYNuDepth)'/Default.VertYNuDepth)./RayCenter((minRow+1):Default.VertYNuDepth,j,3);
                GroundVertEdge(1:(minRow),j) = (1- minRow/Default.VertYNuDepth)./RayCenter(1:(minRow),j,3);
           end
       end
       GroundVertEdgePics = GroundVertEdge;

       TextureFeature.Abs = [ TextureFeature.Abs min([GroundVertEdgePics(:)...
                repmat((Default.VertYNuDepth:-1:1)'/Default.VertYNuDepth,[Default.HoriXNuDepth 1])./reshape(RayCenter(:,:,3),[],1)],[],2)];
%                 /FeaMax(end)];

       TextureFeature.Abs = [SmallSup(:) TextureFeature.Abs];

return;
