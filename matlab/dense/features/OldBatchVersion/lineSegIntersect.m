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
function doIntersect = lineSegIntersect(x,y)

%% we need to find out if the line segment xy(1)xy(2) intersects the line
%% segment xy(3)xy(4)

%% this is true if and only if xy(1) and xy(2) lie on opposite sides of
%% segment xy(3)xy(4) and xy(3) and xy(4) lie on opposite sides of the
%% segment xy(1)xy(2)

%% we use the convex hull trick and check if it is a quadilateral with
%% alternating vertices

disp('lineSegIntersect')
doIntersect=logical(0);
if ((x(2)-x(1))~=0)
    m1 = (y(2)-y(1))/(x(2)-x(1));
    c1 = y(1)-m1*x(1);
    chk1=(((y(3)-m1*x(3)-c1)*(y(4)-m1*x(4)-c1))<0);
else
    chk1=((x(3)-x(1))*(x(4)-x(1))<0);
end

if ((x(4)-x(3))~=0)
    m2 = (y(4)-y(3))/(x(4)-x(3));
    c2 = y(3)-m2*x(3);
    chk2=(((y(1)-m2*x(1)-c2)*(y(2)-m2*x(2)-c2))<0);
else
    chk2=((x(1)-x(3))*(x(2)-x(3))<0);
end

if (chk1 && chk2)
    doIntersect=logical(1);
end

% k=convhull(x,y);
% if length(k)==5
%     if ((k(1)==1)||(k(1)==2))
%         if ((k(2)==3)||(k(2)==4))
%             if ((k(3)==1)||(k(3)==2))
%                 if ((k(4)==3)||(k(4)==4))
%                     doIntersect=logical(1);
%                 end
%             end
%         end
%     elseif ((k(1)==3)||(k(1)==4))
%         if ((k(2)==1)||(k(2)==2))
%             if ((k(3)==3)||(k(3)==4))
%                 if ((k(4)==1)||(k(4)==2))
%                     doIntersect=logical(1);
%                 end
%             end
%         end
%     end
% end
