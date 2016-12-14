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
function [mask,SupIndex,PixelMask,PatchMask]=...
          makeSubSupMask(gridinfo,sup_hi,sup,maskindex,width,hight);

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename;

% calculating the parameter of patch
% [horizontal size of the patch in pixel size, vertical sieze of patch in pixel size]
ratioh(1:2) = gridinfo(:,1)./gridinfo(:,end);
% [horizontal size of the step of point in pixel size, vertical sieze of step of point in pixel size]
ratioh(3:4) = gridinfo(:,1)./gridinfo(:,2);
% [horizontal size of the step of mask in pixel size, vertical sieze of step of mask in pixel size]
MaskStep = ceil(gridinfo(:,2)./gridinfo(:,end));

% initialize the mask
mask = zeros(gridinfo(2,1),gridinfo(1,1)); 
SupIndex = [];
PixelMask = logical(zeros(gridinfo(2,1),gridinfo(1,1)));
PatchMask = logical(zeros(gridinfo(2,2),gridinfo(1,2))); 

for j=(maskindex(1)-MaskStep(1)):MaskStep(1):(gridinfo(1,2)+MaskStep(1)) % j increase in horizontal direction
    for i=(maskindex(2)-MaskStep(2)):MaskStep(2):(gridinfo(2,2)+MaskStep(2)) % i increase in vertical direction
        centerY = ceil(ratioh(4)/2*(2*i-1));
        centerX = ceil(ratioh(3)/2*(2*j-1));
        indicesY = min(max(centerY-hight(2),1),gridinfo(2,1)):max(min(centerY+hight(1),gridinfo(2,1)),1);
        indicesX = min(max(centerX-width(2),1),gridinfo(1,1)):max(min(centerX+width(1),gridinfo(1,1)),1);
        patch_sup = sup_hi( indicesY, indicesX);

        % decide the center subsuperpixel (only do this in the smallest size of patch)
  
        if (j>=1&j<=gridinfo(1,2)&i>=1&i<=gridinfo(2,2))
        	SupTemp = analysesupinpatch(patch_sup,sup(i,j));
            
       		PixelMask(centerY,centerX) = true;
        	PatchMask(i,j) = true;
                SupIndex = [SupIndex; SupTemp];
        else
        	SupTemp = analysesupinpatch(patch_sup);
        end
        mask(indicesY, indicesX) =  patch_sup == SupTemp;
        %SupIndex = [SupIndex; SupTemp];
		
    end
end
