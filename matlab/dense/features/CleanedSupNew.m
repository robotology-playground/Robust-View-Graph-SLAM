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
function [Sup, SupOri, SupNeighborTable]=CleanedSupNew(Default,Sup,maskSky, SupNeighborTable)

displayFlag = false;
ImCloseFlag = false;

% get rid of the Sky for Sup 
SkyCandidate = unique(Sup(maskSky));
% testing sky candidate for two property
% 1) is at least half of the superpixel are sky
tempSkyCandidate = SkyCandidate;
for i = tempSkyCandidate'
    mask = Sup == i;
    if sum(mask(maskSky))/sum(mask(:))<=0.5 % check if half is Sky
       SkyCandidate = setdiff(SkyCandidate', i)';
    end
end

% 2) is at least one of its  neighbor is sky
tempSupNeighborTable = SupNeighborTable;
for i = SkyCandidate' 
    mask = Sup == i;
    Nei = setdiff(find(tempSupNeighborTable(i,:)),i);
    if sum( sum( repmat( SkyCandidate, 1, size(Nei,2)) == ...
                    repmat( Nei, size(SkyCandidate,1), 1) )) >=1% check if at least one neighbor is sky
       Sup(mask) = 0;

       % clearn SupNeighborTable
       SupNeighborTable(i,:) = 0;
       SupNeighborTable(:,i) = 0;
    end
end
maskSky = Sup == 0;

if displayFlag,
DisplaySup(Sup, 2000);
end

SupOri = Sup;
% extend the sky to merger small sup near sky
ThreSmall = 10;
SE = strel('disk',5);
maskSky_dilate = imdilate(maskSky,SE);
SmallSupNearSky = setdiff(Sup(maskSky_dilate),0);
if ~isempty(SmallSupNearSky)
   SE = strel('disk',3);
   for i = SmallSupNearSky'
       mask = Sup ==i;
       if sum(mask(:)) < ThreSmall
          % naive merge
          mask_dilate = imdilate(mask,SE);
          mask_dilate(mask) = 0;
          mask_dilate(maskSky) = 0;
          target = mode(Sup(mask_dilate));
          if isnan(target)
             Sup(mask) = 0;
             % clearn SupNeighborTable
             SupNeighborTable(i,:) = 0;
             SupNeighborTable(:,i) = 0;
          else
             Sup(mask) = target;
             SupNeighborTable( target,:) = ...
                SupNeighborTable( target,:) + SupNeighborTable( i,:);
             SupNeighborTable( :, target) = ...
                SupNeighborTable( :, target) + SupNeighborTable( :, i);
             % clearn SupNeighborTable
             SupNeighborTable(i,:) = 0;
             SupNeighborTable(:,i) = 0;
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


if ImCloseFlag ==1
   % 1) clearn the Sup first (default mask close)
   %    with option to merge under other condition
   %    use just a greedy closing algorithm
   NuSup = setdiff(unique(Sup)',0);
   se = strel('diamond',5);
   for i = NuSup
       mask = Sup == i;
       CloseMask = imclose(mask,se);
       ClosedCandidate = setdiff(unique(Sup(CloseMask)),i);
       if ~isempty(ClosedCandidate)
       for j = ClosedCandidate'
           if all(CloseMask(Sup == j))
              Sup( Sup == j) = i;
           end
       end
       end
   end
end
