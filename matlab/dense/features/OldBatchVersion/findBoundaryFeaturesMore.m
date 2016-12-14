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
function FeaSupNeighList=findBoundaryFeaturesMore(Sup_low,Sup_medi,nList,seglist, FeaMax, slack)
%tic
%load ../debugFeat.mat
global H2;
FeaSize = size(H2,1)

NuSup = setdiff(unique(Sup_low)',0);
NuSup = sort(NuSup);

[currMask1, xVec1, yVec1, currMask2, xVec2, yVec2]=creatingMasks(Sup_low,Sup_medi);
% load /afs/cs/group/reconstruction3d/scratch/maskSup1.mat
% load /afs/cs/group/reconstruction3d/scratch/maskSup2.mat

mask = zeros(size(Sup_medi,1),size(Sup_medi,2));
bMask = mask;
bCurrNeighMask = mask;

%nRows = size(mask,1);
%nCols = size(mask,2);
%newMask = zeros(nRows,nCols,2);
%colCount = 1:nCols;
%rowCount = (1:nRows)';
%newMask(:,:,1) = repmat(colCount,nRows,1);
%newMask(:,:,2) = repmat(rowCount,1,nCols);

currNeighMask=zeros(size(Sup_medi,1),size(Sup_medi,2));
xV=zeros(1,size(mask,2));
yV=zeros(1,size(mask,1));
FeaSupNeighList = zeros(size(NuSup,2),FeaSize);

countSup=1;
for i = NuSup
    nListCurr=find(nList(:,1)==i);
    posInNuSup=find(NuSup==i);
    if (posInNuSup<=500)
        mask = full(currMask1{posInNuSup});
        xV=xVec1{posInNuSup};
        yV=yVec1{posInNuSup};
    else
        mask = full(currMask2{posInNuSup-500});
        xV=xVec2{posInNuSup-500};
        yV=yVec2{posInNuSup-500};
    end
    %% find midpoint of this superpixel
    %% Notice that the sum command sums the columns of a matrix
    x(1)=findSupMid(xV);    
    y(1)=findSupMid(yV);
    
    
    for j=1:length(nListCurr)
        neighNum=nList(nListCurr(j),2);
        posInNuSup=find(NuSup==neighNum);
        if (posInNuSup<=500)
            currNeighMask = full(currMask1{posInNuSup});
             xV=xVec1{posInNuSup};
             yV=yVec1{posInNuSup}; 
        else
            currNeighMask = full(currMask2{posInNuSup-500});
            xV=xVec2{posInNuSup-500};
            yV=yVec2{posInNuSup-500};
        end
        
        %% Now we find the mid-point of the neighbor superpixel
        
        x(2)=findSupMid(xV);
        y(2)=findSupMid(yV);
        
        %iPos=find(FeaSupList(1,:)==i);
        %jPos=find(FeaSupList(1,:)==neighNum);
        
        %Xi = FeaSupList(2:end,iPos);
        %Xj = FeaSupList(2:end,jPos);
        
        %FeaSupNeighList(countSup,3:(3+(length(Xi))-1)) =abs(Xi-Xj)';
        
%         %% Now we calculate the features with just average over the
%         %% boundary
%         [xCen,yCen]=findIntersectionPoint(mask,currNeighMask,x(1:2),y(1:2),newMask);
        xCen = mean(x);
        yCen = mean(y);
        R = [sqrt((xCen-x(1))^2+(yCen-y(1))^2);sqrt((x(2)-xCen)^2+(y(2)-yCen)^2)];
        R=R-slack*R;
        
        bMask = findCloserPointsSquare(mask,xCen,yCen,R(1));
        bCurrNeighMask = findCloserPointsSquare(currNeighMask,xCen,yCen,R(2));
        
%         bMask = findCloserPointsSquare(mask,xCen,yCen,R(1));
%         bCurrNeighMask = findCloserPointsSquare(currNeighMask,xCen,yCen,R(2));
%         
        if ((sum(sum(bMask))==0)||(sum(sum(bCurrNeighMask))==0))
            bXi = zeros(1,FeaSize);
            bXj = zeros(1,FeaSize);
        else
            bXi = (mean(H2(:,logical(bMask))')./FeaMax(1,2:18)); %% Min change this to 18
            bXj = (mean(H2(:,logical(bCurrNeighMask))')./FeaMax(1,2:18)); %% Min change this to 18
        end        
      
        FeaSupNeighList(countSup,1:FeaSize) =abs(bXi-bXj);
        %% Now we check if any of the lines intersects the current
        %% super-pixel and its current neighbor
        FeaSupNeighList(countSup,FeaSize+1)=checkIfALineCrosses(seglist,x,y);
%         for k=1:size(seglist,1)
%             x(3)=seglist(k,1);
%             y(3)=seglist(k,2);
%             x(4)=seglist(k,3);
%             y(4)=seglist(k,4);
%             if(lineSegIntersect(x,y))
%                 bndry_features(countSup,3)=1;
%                 break;
%             end
%         end

%        FeaSupNeighList(countSup,(3+length(Xi)+length(bXi))) = bndry_features(countSup,3);                
%         FeaSupNeighList(countSup,(3+length(Xi))) = bndry_features(countSup,3);                
        countSup=countSup+1;
    end
end

size(FeaSupNeighList)
