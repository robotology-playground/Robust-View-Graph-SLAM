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
function [SupFact, nList] = AnalyzeSup(Sup,maskSky)

% this function analyze the Sup of the NuPatch in each Sup index
NuSup = sort( setdiff( unique(Sup)',0));
edges =[ NuSup];
[yn xn] = size(Sup);

% 1 ) NuPatch in each Sup 
SupFact = []; %NuSup' histc( Sup(:), edges)
nList = [];
for i = NuSup
    mask = Sup == i;
    SupFactTemp = [];
  
    % 2) Center Position of the Sup
    [y x] = find(mask);
    y50 = prctile(y,50);
    x50 = prctile(x,50);
    SupFactTemp = [SupFactTemp y50/yn x50/xn ];

    % 3) x^2 y^2 position of superpixel
    SupFactTemp = [SupFactTemp SupFactTemp(:,(end-1):end).^2 ];

    % 4) x y 10th & 90th
    y90 = prctile(y,90);
    x90 = prctile(x,90);
    y10 = prctile(y,10);
    x10 = prctile(x,10);
    SupFactTemp = [SupFactTemp y90/yn x90/xn y10/yn x10/xn ];

    % 5) eccentricity
    x = x/xn;
    y = y/yn;
    C = cov([x y]);
    [v e] = eig(C);
    tt = diag(e);
    if size(tt,1)~=2
        tt = [tt ;tt];
    end
    [I C] = max(abs(tt));
    ta = v(:,C).*sign(tt(C));
    SupFactTemp = [SupFactTemp sqrt(abs(tt))' acos(ta(1))'];

    % 6) Neighbor count
    SE = strel('disk',3);
    mask_dilate = imdilate(mask,SE);
    mask_dilate_edge = mask_dilate;
    mask_dilate_edge(mask) = 0;
    mask_dilate_edge(maskSky) = 0;
    target =unique(Sup(mask_dilate_edge));
    if any(target <= 0) || any(isnan(target))
       disp('AnalyzeSup error');
    end
    newNei = [i*ones(size(target,1),1) target];
    newNei = sort(newNei,2);
    nList = [ nList;  newNei];
    SupFactTemp = [SupFactTemp size(target,1)];

    SupFact = [SupFact; SupFactTemp];    
end
nList = unique(nList,'rows');
SupFact = [NuSup' histc( Sup(:), edges) SupFact];
return;
