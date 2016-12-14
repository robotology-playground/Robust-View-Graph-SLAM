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
function [GroundBounaryFeature] = GroundBounaryFea(Default, Gmask)

% This function generate the feature related tothe ground boundary position
% From Pinhole camera model, we can derive
% d_g = H * sqrt(f^2+ u_g^2+ v_g^2)./ u_g
% d   = H * sqrt(f^2+ u_g^2+ v_g^2)./ u
% where 
% d_g: depth for ground point
% d : depth for Non ground point
% u_g v_g: image coordinate of poisition of ground point
% u v: image coordinate of poisition of Non ground point
% f: focal length
% H: camera height (to be learned, assuming every image has similiar H)

% Input--
% Gmask: Learned Ground Mask (Used Rajiv Output)

im_x_position = repmat(((1:Default.HoriXNuDepth)-0.5)/Default.HoriXNuDepth - ...
                Default.Ox,[Default.VertYNuDepth 1]);
im_y_position = repmat(((Default.VertYNuDepth:-1:1)'-0.5)/Default.VertYNuDepth - ...
                Default.Oy,[1 Default.HoriXNuDepth]);
RayLength = sqrt(1+ im_x_position.^2+ im_y_position.^2); % Ray length before normalize
im_y_position(~Gmask) = -Inf; % set Non Ground to -Inf to ignore their y_position value
im_y_position( im_y_position >=0) = -Inf; % assume camera center as horizon, so any point higher than horizon just ignore

Ud = zeros(size(Gmask));
for i = 1:Default.VertYNuDepth
    Ud(i,:) = max(im_y_position(i:end,:),[],1);
end

GroundBounaryFeature = abs(RayLength./Ud);

return;
