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
function [FeatureSup ]= f_sup_old(Default,smallSup, sup, SupNeighborTable)
% This fuction calculate the features of each superpixel (total 14)

%%%
% loads a new constant from generalData/SFeaMax.mat
%
% then for each superpixel index, calculate 13 features (each is
% divided by a number in SFeaMax)
% 
% 1) % of image spanned by this superpixel
% 2) x and y position of the center of mass of the superpixel [0-1]
% 3) x and y squared
% 4) x and y positions of 10% and 90% of mass [0-1]
% 5) # of unique superpixels that border this one
% 6) angle of the principal direction of the shape of this
% superpixel (max eigenvector) and sqrt(max eigenvalue)
%%%%

% normalize the Superpixel feature according to the previous result
load([Default.ParaFolder 'SFeaMax.mat']);
% Default.Flag.NormalizeFlag = 1;
imdilateFlag = 0;
% SFeaMax = 10.^floor(log10(SFeaMax));
if  Default.Flag.NormalizeFlag == 1
    SFeaMax = 10.^floor(log10(SFeaMax));
else
    SFeaMax(:) = 1;
end

[yn xn] = size(sup);
NuFea = 13;
NuSup = unique(smallSup(:))';
SupMaxInd = (max(sup(:)));
FeatureSup = sparse(NuFea,SupMaxInd);
SE = strel('diamond',3);
NeighborList = []; % size undetermined
for i=NuSup
%     tic;
    count = 2;
    l = 2;
    FeaturePicsSup = zeros(NuFea,1);
    mask = sup==i;
    % calculating feature
    
    % 2) x y position of superpixel
    [y x] = find(mask);
    y50 = prctile(y,50)/yn;
    x50 = prctile(x,50)/xn;
    FeaturePicsSup(count:(count+1)) = [x50/SFeaMax(1,l); y50/SFeaMax(1,l+1)];
    l = l+2;
    count = count + 2;
    % 3) x^2 y^2 position of superpixel
    FeaturePicsSup(count:(count+1)) = [(x50)^2/SFeaMax(1,l); (y50)^2/SFeaMax(1,l+1)];
    l = l+1;
    count = count + 2;
    % 4) x y 10th & 90th
    y90 = prctile(y,90)/yn;
    x90 = prctile(x,90)/xn;
    y10 = prctile(y,10)/yn;
    x10 = prctile(x,10)/xn;
    FeaturePicsSup(count:(count+3)) = [x10/SFeaMax(1,l); y10/SFeaMax(1,l+1); x90/SFeaMax(1,l+2); y90/SFeaMax(1,l+3)];
    l = l+4;
    count = count + 4;
    % 5) number of connected superpixel
%     if i==56
%        disp('i=56'); 
%     end   
  if imdilateFlag  
    mask_dilate = imdilate(mask,SE);
    mask_dilate_edge = mask_dilate;
    mask_dilate_edge(mask) = 0;
    [list_sup] = unique(sup(mask_dilate_edge)); %hard work
    FeaturePicsSup(count) = size(list_sup,1)/SFeaMax(1,l);

    % storaging the Neighborhood List
    newNei = [i*ones(size(list_sup,1),1) list_sup];
    newNei = sortrows(newNei,2);
    NeighborList = [ NeighborList; newNei];
  end
    lNuNei = l;
    countNuNei =count;
    l=l+1;
    count = count + 1;
    % 6) eccentricity
    [y x] = find(mask);
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
    
    % abs(tt): the standard deviation of the prime axis
    % v(:,C).*sign(tt(C)) : the direction of the prime axis
    FeaturePicsSup(count:(count+2)) =  [sqrt(abs(tt))/SFeaMax(1,l); acos(ta(1))/SFeaMax(1,l+1)];
    if size(FeaturePicsSup,1) ~=13
        disp('error');
    end    
    %l=l+1;
	%MIN: Convert the eigenvector to a positive angle between 0 to 360 
    %FeaturePicsSup = [FeaturePicsSup; sqrt(abs(tt)); abs(v(:,C))]; 
    FeatureSup(:,i) = FeaturePicsSup;
%     toc;
end
Hist = histc( sup(:), NuSup);
FeatureSup(1,NuSup) = Hist/sum(Hist)./SFeaMax(1);
tempNuNeighbor = sum(SupNeighborTable,2);
FeatureSup( countNuNei, NuSup) = tempNuNeighbor(NuSup)./SFeaMax(lNuNei);
