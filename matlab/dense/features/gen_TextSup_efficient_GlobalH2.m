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
function [TextSup]=gen_TextSup_efficient_GlobalH2( param, SelectSegmentationPara)
% process H's superpixels into Hi and Medi Resolution
% this function generate superpixel using default parameter
% but can also change to manually input parameter

% set to global for gen_TextSup_efficient_GlobalH2 only
global H2;

% default parameter
if nargin < 3
    SelectSegmentationPara = 0;
end
DisplayFlag = param.Flag.DisplayFlag; % set to display or not

scale_sigm =[ 1 1.6];
scale_k = [ 1.6 3];
scale_minV = [ 1 3];

%==================== choose 6 different feature channels
 Pick= [1 10 11;
        1 2 5;
        1 3 7;
        10 14 17;
        12 15 13;
        10 10 11];
 NuPick = size(Pick,1);
 reduce = 1; %100 percentage (used to reduce the size to process superpixel)
% ================================

% find the dimension size of the Hi Resolution H
[VertYSizeHiREs,HoriXSizeHiREs,~]= size(H2);
        
% using a median size image to generate superpixel to reduce computation
% intensity (the median size has a upper threshould SegVertYSize SegHoriXSize)
if VertYSizeHiREs*HoriXSizeHiREs > param.SegVertYSize*param.SegHoriXSize

   % Downsample high resolution image to a median size image
   H2 = imresize(H2,([param.SegVertYSize param.SegHoriXSize ]*reduce+4),'nearest'); % +4 because edge error
end
[VertYImg,HoriXImg,~] = size(H2);

%========================================
temp = max(H2,[],1);
temp = max(temp,[],2);
H2 = H2./repmat(temp,[VertYImg,HoriXImg,1]);
%=======================================
    
% Process 6 different feature channel superpixel each with large and median scale
for m=1:NuPick

    img=H2(:,:,Pick(m,:));    
    
    if DisplayFlag
        figure(1); image(img);
    end

    %=================================    
    % choose superpixel of the images
    % default segmentation parameter
    for j = 1:2% number of scale of superpixel
        
        ok = 0; % ok ==1 means accept the segmentation
        while 1
        % call the efficient segment function writen in C++ from MIT
        % Output the high resolution image ( + 1 since the smallest index can be zero) 
        a = segmentImgOpt(param.sigm*scale_sigm(j), param.k*scale_k(j), param.minp*scale_minV(j), uint8(img*255), ...
                           [param.OutPutFolder,param.filename '.ppm'], 0) + 1;
        a = a(2:(end-2),2:(end-2)); % clean the edge superpixel index errors

        % Arrange the superpixel index in order
        %Downsample to size size as prediected depth map
        a = imresize(a,[param.VertYNuDepth param.HoriXNuDepth],'nearest');
        ma = max(a(:));
        Unique_a = unique(a);
        SparseIndex = sparse(ma,1);
        SparseIndex(Unique_a) = 1:size(Unique_a);
        TextSup{m,j} = full(SparseIndex(a));
        clear a SparseIndex Unique_a ma;

        % clean superpixel section ====================================================================
        % merage all small point in higher scale segmentation
%         if j ~= 1
%            TextSup{m,j} = premergAllsuperpixel(TextSup{m,j});
%         end
        % =============================================================================================

        % show superpixel
        if DisplayFlag 
           figure(1);
           imagesc(TextSup{m,j});
           newmap = rand(max(max(TextSup{m,j})),3);
           colormap(newmap);
        end

        % check if need to select segmentation parameter
        if SelectSegmentationPara==1;
           ok = input('Is the segmentation of image OK');% input new segmentation parameter
        else
           ok =1 ;% accept default segmentation parameter
        end

        if ok==1;
           break;
        end

        % Get the user selected parameter
        sigm = input('type sigm of segmentation');
        k = input('type k of segmentation');
        minp = input('type min of segmentation');

        end % end of while 1   
    end % end of j = 1:2 (large and median scale)
end % end of m=1:NuPick (NuPick different feature channel)   

H2 = H2.*repmat(temp,[VertYImg,HoriXImg,1]);
% save([ScratchDataFolder '/data/TextLowResImgIndexSuperpixelSepi' num2str(BatchNu) '.mat'], 'TextLowResImgIndexSuperpixelSep');
return;
