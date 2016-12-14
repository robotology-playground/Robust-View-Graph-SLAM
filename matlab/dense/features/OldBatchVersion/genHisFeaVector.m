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
function [FeaVector] = genHisFeaVector(f,RowTop,RowBottom,ColumnLeft,ColumnRight,i);
% This function set the feature to the right format of Feature Vector
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename NuRow_default;

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

FeaVector = [];
shift = [0 0];
l = 1;
    [Ix Iy] = meshgrid(max(min(ColumnLeft+shift(l,1):ColumnRight+shift(l,1),HoriXNuDepth),1),...
                       max(min(RowTop+shift(l,2):RowBottom+shift(l,2),VertYNuDepth),1));
    maskNeibor = sub2ind([VertYNuDepth, HoriXNuDepth], Iy(:), Ix(:));
    if (sum(sum(Iy==1))== size(Iy(:),1) && l ==4)||(sum(sum( Iy == VertYNuDepth ))== size(Iy(:),1) && l ==5)
        FeaVector =[ FeaVector ;(conv2(f(maskNeibor,:),[0.25; 0.25; 0.25; 0.25],'same'))'];
    else
        FeaVector =[ FeaVector ;f(maskNeibor,:)'];
    end

return;
