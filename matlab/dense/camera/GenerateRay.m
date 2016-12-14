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
function [ray] = GenerateRay(HoriXNuDepth,VertYNuDepth,type,a,b,Ox,Oy);
% This function generate the ray of center or corner of a patch

if nargin < 3
    type = 'center';
    a = 0.70783777 %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
    b = 0.946584169%0.085; % vertical physical size of image plane normalized to focal length (in meter)    Ox = -0.010727086; % camera origin offset from the image center in horizontal direction
    Ox = 0.5; % camera origin offset from the image center in horizontal direction
    Oy = 0.5; % camera origin offset from the image center in vertical direction
elseif nargin < 4
    a = 0.70783777 %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
    b = 0.946584169%0.085; % vertical physical size of image plane normalized to focal length (in meter)    Ox = -0.010727086; % camera origin offset from the image center in horizontal direction
    Ox = 0.5; % camera origin offset from the image center in horizontal direction
    Oy = 0.5; % camera origin offset from the image center in vertical direction
elseif nargin < 5
    b = a;
    Ox = 0.5; % camera origin offset from the image center in horizontal direction
    Oy = 0.5; % camera origin offset from the image center in vertical direction
elseif nargin <6
    Ox = 0.5; % camera origin offset from the image center in horizontal direction
    Oy = 0.5; % camera origin offset from the image center in vertical direction
elseif nargin <7
    Oy = Ox;
end

if type == 'center'
    im_x_position = repmat(((1:HoriXNuDepth)-0.5)/HoriXNuDepth - Ox,[VertYNuDepth 1]);
    im_y_position = repmat(((VertYNuDepth:-1:1)'-0.5)/VertYNuDepth - Oy,[1 HoriXNuDepth]);
elseif type == 'corner'
    im_x_position = repmat((0:(HoriXNuDepth))/(HoriXNuDepth) - Ox,[VertYNuDepth+1 1]);
    im_y_position = repmat((((VertYNuDepth):-1:0)')/(VertYNuDepth) - Oy,[1 HoriXNuDepth+1]);
end

ray = cat(3,a*im_x_position, b*im_y_position, ones(size(im_y_position))); %[ VertYSize horiXSize 3]
ray = ray./repmat(sqrt(sum(ray(:,:,1).^2+ray(:,:,2).^2+ray(:,:,3).^2,3)),[1 1 3]);
return;
