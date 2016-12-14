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
function [AveSupFea ] = AveSupFea(Sup, MedSup, SupFact, MultiScaleSupTable, FeaMax, pick)

% this function generate the AveSupFea for multiScaleSup

global H2;

NuSup = setdiff(unique(Sup)',0);
NuSupSize = size(NuSup,2);
%FeaSup = zeros(NuSupSize,size([H2;H4],1)+1);
%FeaSup = zeros(NuSupSize,size([H2],1)+1);
FeaSup = zeros(NuSupSize,size([H2],1)*3+1); % add up vertical and horizontal features
%size(FeaSup)
l = 1;
for j = NuSup
    mask = MedSup ==j;
    [y x] = find(mask);
    HRange = min(x):max(x);
    VRange = min(y):max(y);
    if sum(mask(:))==0
       disp(' AveSupFea error');
    end
%    FeaSup(l,:) = [j (mean(H2(:,mask)')./FeaMax(1,2:18)) (mean(H4(:,mask)')./FeaMax(1,19:35))];
    if pick == 1
       FeaSup(l,:) = [j (mean(H2(:,mask)')./sqrt(FeaMax(1,2:18))) ...
                        (mean( reshape( H2(:,VRange,:), size(H2,1), [])')./sqrt(FeaMax(1,2:18))) ...
                        (mean( reshape( H2(:,:,HRange), size(H2,1), [])')./sqrt(FeaMax(1,2:18)))]; % keeping the Sup index at the first column in FeaSup
    elseif pick == 2
       FeaSup(l,:) = [j (mean(H2(:,mask)')./FeaMax(1,2:18)) ...
                        (mean( reshape( H2(:,VRange,:), size(H2,1), [])')./(FeaMax(1,2:18))) ...
                        (mean( reshape( H2(:,:,HRange), size(H2,1), [])')./(FeaMax(1,2:18)))]; % keeping the Sup index at the first column in FeaSup
    else
       FeaSup(l,:) = [j (mean(H2(:,mask)')./FeaMax(1,19:35)) ...
                        (mean( reshape( H2(:,VRange,:), size(H2,1), [])')./(FeaMax(1,19:35))) ...
                        (mean( reshape( H2(:,:,HRange), size(H2,1), [])')./(FeaMax(1,19:35)))]; % keeping the Sup index at the first column in FeaSup
    end   
    l = l + 1;
end

disp('go to AveSupFea');
% calculate the MultiScale Fea of Sup
FeaSize = (size(H2,1));
ScaleSize = size(MultiScaleSupTable,2)-1;
AveSupFea = zeros(NuSupSize,FeaSize*ScaleSize);
l = 1;
for j = NuSup
    row = find(MultiScaleSupTable(:,1) == j);
    Target = MultiScaleSupTable( row, 2:end);
    MultiMask = (MultiScaleSupTable(:,2:end) == repmat(Target,[ NuSupSize 1]));
    ScaleMask = sum(MultiMask,1) > 1; % if Scale have more than one Sup in the save Scale
    MultiMask(row,ScaleMask) = false; % set the j to zeros
    Wei = repmat(SupFact(:,2),[1 ScaleSize]).*MultiMask; % get rid of the j sup
    Wei = Wei ./ repmat(sum(Wei,1),[NuSupSize 1]); % normalize the Wei to have sum to 1
    AveSupMultScaleFea = sum(repmat(FeaSup(:,2:(FeaSize+1)),[1 1 ScaleSize]).*repmat(permute(Wei,[1 3 2]),[1 FeaSize 1]),1);
    AveSupFea(l,:) = AveSupMultScaleFea(:)';
    l = l + 1;
end
AveSupFea = [FeaSup AveSupFea];

return;


