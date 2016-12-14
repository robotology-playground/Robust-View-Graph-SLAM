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
function [f, fInd] = AbsFeatureGenMex_MemoryEfficient(param, SmallSup, HiSupi, SupMaskFlag, FeaMax, fInd)

% This function generate the average of the feature within a certain mask of a sample point
% Input--
% H2: texture filter output
% HiSup: the Hi resolution superpixel index matrix
% SupMaskFlag: if SupMaskFlag is 1, we calculate the averaging using superpixel mask (irregular mask)
% Output--
% f: Feature matrix of size(No od depth point, No of feature vector)

global H2;

% Parameter
[ImgResY,ImgResX] = size(H2(:,:,1));
   
% define image and sample point and patch size infomation
gridinfo = [param.TrainHoriXSize param.HoriXNuDepth param.HoriXNuPatch; param.TrainVerYSize param.VertYNuDepth param.VertYNuPatch];

% Grid Info
ratio(1:2) = floor(gridinfo(:,1)./gridinfo(:,end) );
ratio(3:4) = floor( gridinfo(:,1)./gridinfo(:,2) );
   
% Patch shape infomation
%hcol = ones(floor(ratio(2)),1);
%hrow = ones(1,floor(ratio(1)));

if SupMaskFlag == 1 % Need to use irregular Superpixel Mask

    % calculate how many mask we need
   NuMask = ceil(gridinfo(:,2)./gridinfo(:,3));
   
   % calcuate the position of the mask
   hight(1) = round((ratio(2)-1)/2);
   hight(2) = ratio(2) - 1 - hight(1);
   width(1) = round((ratio(1)-1)/2);
   width(2) = ratio(1) - 1 - width(1); 

   row_start = 1;
   f = [];
   f_pics_mask = zeros( gridinfo(1,3)*gridinfo(2,3), 17);
   for j = 1:NuMask(1)
       for k = 1:NuMask(2);
           fIndNew = fInd;
           % first generate the mask respect to the dominate subsuperpixel
           [mask,PixelMask,PatchMask] = makeSubSupMaskNew2(gridinfo, HiSupi, SmallSup, [j; k], width, hight);
%           [mask,SupIndex,PixelM,PatchMask] = makeSubSupMask(gridinfo, HiSupi, SmallSup, [j; k], width, hight);
            

           % calculaing the normalize value
 %          NormalizeValue = conv2(hcol,hrow,mask,'same');%///////////////////////////////////////////////

           % generate the 1:34 features for H2  for 1 center and 4 neighbor (left right top bottom)
           for m = 1:17
               
%               temp = conv2(hcol, hrow, H2(:,:,m).*mask, 'same');%///////////////////////////////////////////
%               tt = temp(PixelM)./NormalizeValue(PixelM);
%               vv =SparseAverageSample2D(H2(:,:,m), floor(ratio(2)), floor(ratio(1)),PixelMask,mask);
%               [floor(ratio(1)), floor(ratio(2))]
               f_pics_mask(:,m) = ...
                              SparseAverageSample2DOptimized(H2(:,:,m),...
                              ratio(2), ratio(1), ...
                              PixelMask,double(mask))...
                             ./FeaMax(1,fIndNew);
               fIndNew = fIndNew+1;
           end

           f(PatchMask, row_start:row_start+size(f_pics_mask,2)-1) = f_pics_mask;
       end
   end
   fInd = fIndNew;

else
   DepthGridSizeY = ImgResY/param.VertYNuDepth;
   DepthGridSizeX= ImgResX/param.HoriXNuDepth;
   % 1) generating the PixelMask
%   PixelMask = logical(zeros(ImgResY,ImgResX));
   [X,Y] = meshgrid(ceil((1/2)*DepthGridSizeX:DepthGridSizeX:ImgResX),...
                    ceil((1/2)*DepthGridSizeY:DepthGridSizeY:ImgResY));
   %PixelMask = sub2ind(size(PixelMask),Y(:),X(:));
   for m = 1:17
       f(:,m) = SparseAverageSample2DOptimized(H2(:,:,m),ratio(2),ratio(1),[Y(:) X(:)], double(ones(size(H2(:,:,m)))))...
                        ./FeaMax(1,fInd);
       fInd = fInd +1;
   end

end

return; 
