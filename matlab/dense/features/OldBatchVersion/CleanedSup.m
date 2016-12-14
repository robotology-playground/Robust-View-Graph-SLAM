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
function [Sup,MedSup, status]=CleanedSup(Sup,MedSup,maskSky)

displayFlag = false;
status = 0;

% This function clean the sup up to certain cretiria
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName;

% 0)get ride of all the nonconnedted superpixel in MedSup
NuMedSup = unique(MedSup)';
SE = strel('disk',3);
for i = NuMedSup
    mask = MedSup ==i;
    [L,num] = bwlabel(mask,8);
    n = histc(L(:),1:num);
    [V maxL] = max(n);
    for j = setdiff(1:num,maxL)
        mask_dilate = imdilate(L == j,SE);
        mask_dilate(mask) = 0;
        if all(mask_dilate(:) == 0)
           %j
           disp('error using Mode empty')
        end
        MedSup(L ==j) = mode(double(MedSup(mask_dilate)));
    end
end

% 1) Make the MedSup have the same index as Sup
NuSup = unique(Sup)';
NuMedSup = unique(MedSup)';
NuResidualSup = setdiff(NuMedSup,NuSup);
if ~isempty(NuResidualSup)
    [Sup, MedSup]=Medi2SupIndex(Sup,MedSup);
end
%check
NuSup = unique(Sup)';
NuMedSup = unique(MedSup)';
NuResidualSup = setdiff(NuMedSup,NuSup);
if ~isempty(NuResidualSup)
   disp('error in CleanSup');
end

% get rid of the Sky for both Sup and MedSup
SkyCandidate = unique(Sup(maskSky));
% check connecting for SkyResidual
% Do we need to deal with these residual??????
for i = SkyCandidate' 
    mask = Sup == i;
    if sum(mask(maskSky))/sum(mask(:))>0.5
       Medmask = MedSup == i;
       Sup(mask) = 0;
       MedSup(Medmask) = 0; 
    end
end
maskSky = Sup == 0;

if displayFlag,
DisplaySup(Sup, 2000);
end

% extend the sky to merger small sup near sky
ThreSmall = 10;
SE = strel('disk',5);
maskSky_dilate = imdilate(maskSky,SE);
SmallSupNearSky = setdiff(Sup(maskSky_dilate),0);
if ~isempty(SmallSupNearSky)
for i = SmallSupNearSky'
    mask = Sup ==i;
    if sum(mask(:)) < ThreSmall
       Medmask = MedSup == i;
       % naive merge
       SE = strel('disk',3);
       mask_dilate = imdilate(mask,SE);
       mask_dilate(mask) = 0;
       mask_dilate(maskSky) = 0;
       target = mode(Sup(mask_dilate));
       if isnan(target)
          Sup(mask) = 0;
          MedSup(Medmask) = 0;
       else
          Sup(mask) = target;
          MedSup(Medmask) = target;
       end
    end
end
end

if displayFlag,
figure(250), imagesc(Sup),
newmap = [rand(max(Sup(:)),3); [0 0 0]];
colormap(newmap);
%DisplaySup(Sup, 3500);
end

% 1) clearn the Sup first (default mask close)
%    with option to merge under other condition
%    use just a greedy closing algorithm
NuSup = setdiff(unique(Sup)',0);
for i = NuSup
    Medmask = MedSup == i;
    se = strel('diamond',10*3);
    CloseMask = imclose(Medmask,se);
    ClosedCandidate = setdiff(unique(MedSup(CloseMask)),i);
    if ~isempty(ClosedCandidate)
    for j = ClosedCandidate'
        if all(CloseMask(MedSup == j))
           Sup( Sup == j) = i;
           MedSup( MedSup == j) = i;
        end
    end
    end
end

%Final check
NuSup = unique(Sup)';
NuMedSup = unique(MedSup)';
NuResidualSup = setdiff(NuMedSup,NuSup);
if ~isempty(NuResidualSup)
   disp('error in CleanSup');
%   status = 1;
end

