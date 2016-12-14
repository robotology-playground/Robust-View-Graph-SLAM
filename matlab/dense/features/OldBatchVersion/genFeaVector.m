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
function [FeaVector] = genFeaVector(f,fsup,VList,HList,i,near,NeighborFeaList);
%function [FeaVector] = genFeaVector(f,fsup,RowTop,RowBottom,ColumnLeft,ColumnRight,i,NeighborFeaList);
% This function set the feature to the right format of Feature Vector
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename NuRow_default;

if nargin < 7 
	NeighborFeaList = 1:5;
end

NuRow = NuRow_default;

% load pics info
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

% generate the range of the row for the same thi (weight value)
%RowskyBottom = floor(NuRow/2);
%PatchSkyBottom = ceil(VertYNuDepth*(1-Horizon));
%if row <= RowskyBottom
%   PatchRowRatio = PatchSkyBottom/RowskyBottom;
%   RowTop = round((row-1)*PatchRowRatio+1);
%   RowBottom = round(row*PatchRowRatio);
%else
%   PatchRowRatio = (VertYNuDepth-PatchSkyBottom)/(NuRow-RowskyBottom);
%   RowTop = round((row-RowskyBottom-1)*PatchRowRatio+1)+PatchSkyBottom;
%   RowBottom = round((row-RowskyBottom)*PatchRowRatio)+PatchSkyBottom; 
%end

FeaVector = [];
% Superpixel Feature
if size(fsup,2)~=0
   shift = [0 0; -1 0; 1 0; 0 -1; 0 1]; % left right top bottom
   for l = 1:1
      [Ix Iy] = meshgrid(max(min(HList+shift(l,1),HoriXNuDepth),1),...
                     max(min(VList+shift(l,2),VertYNuDepth),1));
      maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
      if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4) || (sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
         FeaVector =[ FeaVector ;(conv2(fsup(:,f(maskNeibor,1)),[0.25; 0.25; 0.25; 0.25],'same'))];
      else
         FeaVector =[ FeaVector ;fsup(:,f(maskNeibor,1))];
      end
   end
end
% Hi Resolution
shift = [0 0; -1 0; 1 0; 0 -1; 0 1];

for l = NeighborFeaList
    [Ix Iy] = meshgrid(max(min(HList+shift(l,1),HoriXNuDepth),1),...
                       max(min(VList+shift(l,2),VertYNuDepth),1));
    maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
    if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4)||(sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
        FeaVector =[ FeaVector ;(conv2(f(maskNeibor,2:35),[0.25; 0.25; 0.25; 0.25],'same'))'];
    else
        FeaVector =[ FeaVector ;f(maskNeibor,2:35)'];
    end
end

% Med 1/3 Resolution
if near ~=1
   shift = 3*[0 0; -1 0; 1 0; 0 -1; 0 1];
end
for l = NeighborFeaList
    [Ix Iy] = meshgrid(max(min(HList+shift(l,1),HoriXNuDepth),1),...
                       max(min(VList+shift(l,2),VertYNuDepth),1));
    maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
    if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4)||(sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
        FeaVector =[ FeaVector ;(conv2(f(maskNeibor,36:69),[0.25; 0.25; 0.25; 0.25],'same'))'];
    else
        FeaVector =[ FeaVector ;f(maskNeibor,36:69)'];
    end
end

% Low 1/9 REsolution
if near ~=1
   shift = 9*[0 0; -1 0; 1 0; 0 -1; 0 1];
end
for l = NeighborFeaList
    [Ix Iy] = meshgrid(max(min(HList+shift(l,1),HoriXNuDepth),1),...
                       max(min(VList+shift(l,2),VertYNuDepth),1));
    maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
    if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4)||(sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
        FeaVector =[ FeaVector ;(conv2(f(maskNeibor,70:103),[0.25; 0.25; 0.25; 0.25],'same'))'];
    else
        FeaVector =[ FeaVector ;f(maskNeibor,70:103)'];
    end
end

% other features
shift = [0 0];%; -1 0; 1 0; 0 -1; 0 1];
for l = 1:1
    [Ix Iy] = meshgrid(max(min(HList+shift(l,1),HoriXNuDepth),1),...
                       max(min(VList+shift(l,2),VertYNuDepth),1));
    maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
    %if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4)||(sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
    %    maskNeibor2 = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:)+1, Ix(:));
    %    FeaVector =[ FeaVector ;(f(maskNeibor,104)+f(maskNeibor2,104))'/2];
    %else
        FeaVector =[ FeaVector ;f(maskNeibor,104)'];
    %end
end

