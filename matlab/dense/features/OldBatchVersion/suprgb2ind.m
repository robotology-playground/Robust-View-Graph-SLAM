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
function [imind, list] = suprgb2ind(imrgb);
% This function take superpixel images in rgb form then output the superpixel
% images in index image form - a yn * xn matrix with superpixel
% indicies as values

% parameter
[yn xn t]= size(imrgb);%find the dimension size of the data matrix
list =[];
imrgb = permute(imrgb,[3 1 2]);
imrgb = imrgb(:,:); % imrgb is now 3 x (yn*xn)
unind = logical(ones(1,yn*xn));
ind = 0;
while any(unind)
    temp = imrgb(:,unind);
    tt = temp(:,1); %first unindexed pixel
    % update unind, d=vector same size as unind w/ 1 if pixel==tt,
    % 0 o/w
    d = (temp(1,:)==tt(1)&temp(2,:)==tt(2)&temp(3,:)==tt(3));
    tar = logical(zeros(1,yn*xn));
    % tar is length of image, w/ 1s from d
    tar(:,unind) = d;
    unind(:,unind)=(~d);
    % make list
    ind = ind + 1;
    list = [ list [ind;double(tt);sum(d)]];
    imind(tar) = ind;
end            
imind = reshape(imind,[yn xn]);
return;

