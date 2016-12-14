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
%function []=CalFeatures(FeaBatchNumber,HistFeaType,AbsFeaType);
function []=CalFeatures(taskName,ImgFolder,TrainSet,LearnType,LearnSkyEx,LearnLog,LearnNear,...
                   LearnAlg,LearnDate,AbsFeaType,AbsFeaDate,HistFeaType,HistFeaDate, ...
                   GeneralDataFolder, ScratchDataFolder, LocalFolder, ClusterExecutionDirectory, FeaBatchNumber)
% This funciton is the meta funciton to choose different feature calculation method
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName...
    TrainVerYSize TrainHoriXSize MempryFactor;




%HistFeaType
%AbsFeaType
%switch AbsFeaType
%   case 'Whole'
%gen_feature_filt1D_sepH2H4_fixMem(FeaBatchNumber,HistFeaType,1);
%      disp('AbsFeaType is Whole');
gen_feature_general(taskName,ImgFolder,TrainSet,LearnType,LearnSkyEx,LearnLog,LearnNear,...
                   LearnAlg,LearnDate,AbsFeaType,AbsFeaDate,HistFeaType,HistFeaDate, ...
                   GeneralDataFolder, ScratchDataFolder, LocalFolder, ClusterExecutionDirectory, FeaBatchNumber,1);
   %case 'Sub'
   %   disp('AbsFeaType is Subsuperpixel');
   %    gen_feature_sep(FeaBatchNumber,HistFeaType,1);      
%   case 'WholeH13'	%Calculate Histogram based features
%      disp('AbsFeaType is WholeH13');
%      gen_feature_H1_H3(FeaBatchNumber,HistFeaType,1);
   %otherwise
   %   disp('AbsFeaType is None.');
   %    gen_feature_sep(FeaBatchNumber,HistFeaType,0);
%end
return;
