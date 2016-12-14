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
function [FeatureSup]= f_sup(sup)
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

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch HoriXNuPatch a_default b_default Ox_default Oy_default...
    Horizon_default filename;

% normalize the Superpixel feature according to the previous result
load([GeneralDataFolder '/SFeaMax.mat']);
SFeaMax = 10.^floor(log10(SFeaMax));

FeatureSup = [];
[yn xn] = size(sup);
NuSup = max(max(sup));
for i=1:NuSup
%     tic;
    l = 1;
    FeaturePicsSup = [];
    mask = sup==i;
    % calculating feature
    % 1) size of superpixel normalized to the total number of subsuperpixel
    NuSub = sum(sum(mask));
    FeaturePicsSup = [FeaturePicsSup; NuSub/(xn*yn)/SFeaMax(1,l)];
    l=l+1;
    % 2) x y position of superpixel
    NuSubSupVertY =  sum(mask,1);
    PosiHoriX = find(NuSubSupVertY);
    AccuHoriX =cumsum(NuSubSupVertY(1,NuSubSupVertY~=0)); 
    x = PosiHoriX(min(find(AccuHoriX>=round(NuSub/2))))/xn;% normalize by the size of the pics
    NuSubSupHoriX =  sum(mask,2);
    PosiVertY = find(NuSubSupHoriX);
    AccuVertY =cumsum(NuSubSupHoriX(NuSubSupHoriX~=0,1)); 
    y = PosiVertY(min(find(AccuVertY>=round(NuSub/2))))/yn;% normalize by the size of the pics
    FeaturePicsSup = [FeaturePicsSup; x/SFeaMax(1,l); y/SFeaMax(1,l+1)];
    l = l+2;
    % 3) x^2 y^2 position of superpixel
    FeaturePicsSup = [FeaturePicsSup; x^2/SFeaMax(1,l); y^2/SFeaMax(1,l+1)];
    l = l+1;
    % 4) x y 10th & 90th
    xten = PosiHoriX(min(find(AccuHoriX>=round(NuSub/10))))/xn;% normalize by the size of the pics    
    xninety = PosiHoriX(min(find(AccuHoriX>=round(NuSub*9/10))))/xn;% normalize by the size of the pics
    yten = PosiVertY(min(find(AccuVertY>=round(NuSub/10))))/yn;% normalize by the size of the pics
    yninety = PosiVertY(min(find(AccuVertY>=round(NuSub*9/10))))/yn;% normalize by the size of the pics
    FeaturePicsSup = [FeaturePicsSup; xten/SFeaMax(1,l); yten/SFeaMax(1,l+1); xninety/SFeaMax(1,l+2); yninety/SFeaMax(1,l+3)];
    l = l+4;
    % 5) number of connected superpixel
%     if i==56
%        disp('i=56'); 
%     end    
    SE = strel('diamond',3);
    mask_dilate = imdilate(mask,SE);
    mask_dilate_edge = mask_dilate;
    mask_dilate_edge(mask) = 0;
    [list_sup] = unique(sup(mask_dilate_edge)); %hard work
    FeaturePicsSup = [FeaturePicsSup; size(list_sup,1)/SFeaMax(1,l)];
    l=l+1;
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
    FeaturePicsSup = [FeaturePicsSup; sqrt(abs(tt))/SFeaMax(1,l); acos(ta(1))/SFeaMax(1,l+1)];
    if size(FeaturePicsSup,1) ~=13
        disp('error');
    end    
    %l=l+1;
	%MIN: Convert the eigenvector to a positive angle between 0 to 360 
    %FeaturePicsSup = [FeaturePicsSup; sqrt(abs(tt)); abs(v(:,C))]; 
    FeatureSup = [FeatureSup FeaturePicsSup];
%     toc;
end
