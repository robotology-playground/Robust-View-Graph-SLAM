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
function [] = RunCleanedSup(BatchNu)

% this function generate the CleanedSup in /data/CleanedSup directory

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName;

load([ScratchDataFolder '/data/LowResImgIndexSuperpixelSep.mat']);
load([ScratchDataFolder '/data/MaskGSky.mat']);

BatchSize = 10;
NuPics = size(filename,2);
BatchRow = 1:BatchSize:NuPics;
% ImageVector = 1:NuPics
%ImageVector = [11 12 14 15 21 22 24 47 48 57 59 62 64 65 70 71 73 74 76 77 83 84 88 90 96 97 99 100 107 112 115 119 120 124 129 132 133 136 138 140 141 145 146 147 154 155  157 159 164 165 168 169 176 178 182 190 194 195 199 202 206 207 208 209 210 212 214 216 217 222 224 320] % Pics form occlusion learning
for i = BatchRow(BatchNu):min(BatchRow(BatchNu)+BatchSize-1,NuPics)
%for i = 1:NuPics
    %i
    load([ScratchDataFolder '/data/MedSeg/MediResImgIndexSuperpixelSep' num2str(i) '.mat']);
    [Sup,MedSup , status]=CleanedSup(LowResImgIndexSuperpixelSep{i},MediResImgIndexSuperpixelSep,maskSky{i});
    if status
       return;
    end
    save([ScratchDataFolder '/data/CleanSup/CleanSup' num2str(i) '.mat'],'Sup','MedSup');
end
