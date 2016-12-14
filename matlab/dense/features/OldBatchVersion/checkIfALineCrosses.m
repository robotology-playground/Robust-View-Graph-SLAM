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
function bndryFeat=checkIfALineCrosses(seglist,x,y)

% This function called by findBoundaryFeaturesMore
% written by Rajiv

x3=seglist(:,1);
y3=seglist(:,2);
x4=seglist(:,3);
y4=seglist(:,4);

bndryFeat=logical(0);

if ((x(2)-x(1))~=0)
    m1 = (y(2)-y(1))/(x(2)-x(1));
    c1 = y(1)-m1*x(1);
    chk1=(((y3-m1*x3-c1).*(y4-m1*x4-c1))<0);
else
    chk1=((x3-x(1)).*(x4-x(1))<0);
end

vertChk=abs((x4-x3));
vertLines=find(vertChk==0);
nonVertLines=find(vertChk~=0);

x3v=seglist(vertLines,1);
y3v=seglist(vertLines,2);
x4v=seglist(vertLines,3);
y4v=seglist(vertLines,4);

x3nv=seglist(nonVertLines,1);
y3nv=seglist(nonVertLines,2);
x4nv=seglist(nonVertLines,3);
y4nv=seglist(nonVertLines,4);

m2 = (y4nv-y3nv)./(x4nv-x3nv);
c2 = y3nv-m2.*x3nv;
chk2(nonVertLines,1)=(((y(1)-m2*x(1)-c2).*(y(2)-m2*x(2)-c2))<0);
chk2(vertLines,1)=((x(1)-x3v).*(x(2)-x3v)<0);

if (sum(chk1 & chk2))
    bndryFeat=logical(1);
end
