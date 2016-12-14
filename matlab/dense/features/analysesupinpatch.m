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
function [maxindex] = analysesupinpatch(patch_sup,GuessIndex);

%%%
% looks like it just outputs the index that is most common in
% patch_sup 
%%%

% Confirmed it just outputs the index that is most common in
% By Min Dec 20 2006
if nargin < 2
    GuessIndex = 0;
end
patch_sup = patch_sup(:);
list = [];
mask = logical(ones(size(patch_sup)));

% try the Guess first
nPixels = sum(patch_sup == GuessIndex);
if nPixels > length(patch_sup)/2
	maxindex = GuessIndex;
	return;
end

% Search All
nVector = unique(patch_sup)'; %MIN: Use indices which are in low-res image.
for n = nVector
	nPixels = sum(patch_sup == n);
	list = [list, nPixels];
	if nPixels > length(patch_sup)/2
		break;
        end
end
[V ind] = max(list);
maxindex = nVector(ind);

return;

