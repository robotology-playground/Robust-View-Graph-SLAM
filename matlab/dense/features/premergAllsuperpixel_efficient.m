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
function [im, SupNeighborTable]=premergAllsuperpixel_efficient(im, Default)
%(This is program conver a non_ordered non_connected superpixel image to 
% ordered superpixel image
% input:
% im = N by M matrix depends on the image size       
% Default - Default.SmallThre = the samllest Sup size

% output:
% im_order = N by M matrix which is ordered

%%%
% For each of the superpixel indicies, finds all of the
% disconnected components with that label.  if the component is
% small (< 200 pixels) or isn't the biggest component with that
% index, call analysesupinpatch with the outline of the component
% to replace it's label with the one that is most common in the outline
%%%

if nargin <2
   Default.SmallThre = 5; %smallest sup size
end
SupNeighborTableFlag = 1;

[yn xn] = size(im);
NuSup = unique(im(:))';
SE = strel('octagon',3);  
for i=NuSup

        % label connected component
        temp = zeros(size(im));
        temp(im(:,:)==i)=1;
        [L,num] = bwlabel(temp,4);

        % find the main piece
         [maxL dum]= mode(L(L~=0));
%        his = histc(L(:), 1:num);
%        [dum maxL ]= max(his);

           if dum > Default.SmallThre;
              SupMerg = setdiff(1:num,maxL);
           else
              SupMerg = 1:num;
           end

           for k = SupMerg
               mask = L==k;
               % then assign those pixels to mostlikely 3 by 3 neighborhood
                mask_dilate = imdilate(mask,SE);
%                mask_dilate = mask | [zeros(yn,1) mask(:,1:(end-1))] ...
%                                   | [mask(:,2:(end)) zeros(yn,1)] ...
%                                   | [zeros(1,xn) ;mask(1:(end-1),:)] ...
%                                   | [mask(2:(end),:); zeros(1, xn)] ....
%                                   | [[zeros(yn-1,1) mask(2:end,1:(end-1))]; zeros(1,xn)]...
%                                   | [[mask(2:end,2:(end)) zeros(yn-1,1) ]; zeros(1,xn)]...
%                                   | [zeros(1,xn); [zeros(yn-1,1) mask(1:(end-1),1:(end-1))]]...
%                                   | [zeros(1,xn); [mask(1:(end-1),2:(end)) zeros(yn-1,1)]]...
%                                   ;                          
               mask_dilate(mask) = 0;
%                im(mask) = analysesupinpatch(im(mask_dilate));%hard work
               im(mask) = mode(im(mask_dilate));
           end
end

% merge the small superpixel with the surrrounding one if it's neighbor is only one
 MaxSupIndex = max(NuSup(:));
 SupNeighborTable = sparse(MaxSupIndex,MaxSupIndex);

if SupNeighborTableFlag
   for i = 1:((xn-1)*yn)
     
       % registed the neoghbor in right
       SupNeighborTable(im(i),im(i+yn)) = 1;
       SupNeighborTable(im(i+yn),im(i)) = 1;

       % registed the neoghbor in below
       if mod(i,yn) ~=0
          SupNeighborTable(im(i),im(i+1)) = 1;
          SupNeighborTable(im(i+1),im(i)) = 1;
       end
   end

   % find out the single neighbor ones and merge them with neighbors
   SingleTar = sum( SupNeighborTable,1);
   for i = find(SingleTar == 1)
       mask = im == i;
       im(mask) = find(SupNeighborTable(:,i) == 1);
       SupNeighborTable(:,i) = 0;
       SupNeighborTable(i,:) = 0;
   end

end
return;
