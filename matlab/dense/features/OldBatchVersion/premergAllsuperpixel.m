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
function [im]=premergsuperpixel(im)
%(This is program conver a non_ordered non_connected superpixel image to 
% ordered superpixel image
% input:
% im = N by M matrix depends on the image size       
% Position3DTure =3 by Q matrix the 3 column entries are [x y z]' in 3d coordinate
% output:
% im_order = N by M matrix which is ordered

%%%
% For each of the superpixel indicies, finds all of the
% disconnected components with that label.  if the component is
% small (< 200 pixels) or isn't the biggest component with that
% index, call analysesupinpatch with the outline of the component
% to replace it's label with the one that is most common in the outline
%%%

if nargin < 2
    rankcheck = 0;
end
[yn xn] = size(im);
Ngmax = max(max(im)); % number of superpixel group
SE = strel('octagon',3);  
for i=1:double(Ngmax)
    if sum(sum(im(:,:)==i))
        % label connected component
        temp = zeros(size(im));
        temp(im(:,:)==i)=1;
        [L,num] = bwlabel(temp,4);
        % find the main piece
        for k=1:num
            his(1,k)=sum(sum(L==k));
        end
        [dum maxL ]= max(his);
        for k=1:num
            mask = L(:,:)==k;
            % filter out the superpixel smaller than 5
            % then assign those pixels to mostlikely 3 by 3 neighborhood
            if k~=maxL || sum(sum(L(:,:)==k))<200
                    mask_dilate = imdilate(mask,SE);
                    mask_dilate(mask) = 0;
                    im(mask) = analysesupinpatch(im(mask_dilate));%hard work
%                     [list_sup] = analysesupinpatch(im(mask_dilate));%hard work
%                     [max_num max_sup_index] = max(list_sup(2,:));
%                     im(mask)=list_sup(1,max_sup_index);
            else
                im(mask) = i;
            end    
        end    
    end    
end
return;
