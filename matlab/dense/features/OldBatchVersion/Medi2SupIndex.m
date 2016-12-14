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
function [Sup, MedSup]=Medi2SupIndex(Sup,MedSup);

% This funciton Merge the superpixel in MediResImgIndexSuperpixelSep to
% have the same unique index of LowResImgIndexSuperpixelSep 
% We need to do this because LowResImgIndexSuperpixelSep is the down-sample of 
% MediResImgIndexSuperpixelSep()

% define global variable
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize;

% initalize parameters
NuSup = unique(Sup)';
NuMedSup = unique(MedSup)';
NuResidualSup = setdiff(NuMedSup,NuSup);

maskResidual = zeros(size(MedSup));
for i = NuResidualSup
    disp(i);
%    size(maskResidual)
%    size(MedSup)
    maskResidual = maskResidual | MedSup ==i;
end

% Start checking if rand of Sup >= RankThresh
for i = NuResidualSup
    mask = MedSup==i;
    SE = strel('diamond',3);
    mask_dilate = imdilate(mask,SE);
    mask_dilate_edge = mask_dilate;
    mask_dilate_edge(maskResidual) = 0;
    MergeTarget = mode(double(MedSup(mask_dilate_edge)));
    if ~any(MergeTarget == NuSup)
        disp('error');
        return
    end
    MedSup(mask) = MergeTarget;
end

return;
