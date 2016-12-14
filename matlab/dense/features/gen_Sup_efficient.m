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
function [MedSup, Sup, Default SupNeighborTableFlag]=gen_Sup_efficient(Default, img, SelectSegmentationPara);
% this function generate superpixel using default parameter
% but can also change to manually input parameter

%%% Jeff's Comments (Min modified)
% Send that image to the CMU
% segmentation program with params 0.8*[sigm, k, min].  If
% SelectSegmenationPara is true, then display the results and ask
% the user to enter new sigm, k, and min; repeat until user is happy.
%
% Output the MedSup(only the smallest scale) and Sup(three scale)
%%%%

% default parameter
if nargin < 3
    SelectSegmentationPara = 0; % if SelectSegmentationPara == 1, enable the parameter interation with user.
end
DisplayFlag = Default.Flag.DisplayFlag; % set to display or not

scale =[0.8 1.6 5]; % use different scale to generate small(0.8) middle(1.6) 5(large) scale of superpixel

    [VertYSizeHiREs,HoriXSizeHiREs,dummy]= size(img);% find the dimension size of the Hi Resolution image
    clear dummy;

    % using a fixed range of median size image [SegVertYSize SegHoriXSize ] 
    %  to generate superpixel to reduce computation
    if VertYSizeHiREs*HoriXSizeHiREs > Default.SegVertYSize*Default.SegHoriXSize

       % Downsample high resolution image to a fixed median size image
       %**** +4 because segmentImgOpt gives constant additinal rows and column
       %**** so add 4 rows and columns a prior then delete then at line 55 
       img = imresize(img,[Default.SegVertYSize+4 Default.SegHoriXSize+4 ],'nearest'); 
    else
       Default.SegVertYSize = VertYSizeHiREs-4;
       Default.SegHoriXSize = HoriXSizeHiREs-4;
    end

% generate superpixel of each image
for j = 1:3% number of scale of superpixel
            
    % choose superpixel of the images
    % default segmentation parameter
    ok = 0; % ok ==1 means accept the segmentation
    while 1

        % call the efficient segment function writen in C++ from MIT
        % Output the high resolution image ( + 1 since the smallest index can be zero)
        if j ==1
	   a = segmentImgOpt( Default.sigm*scale(j), Default.k*scale(j), Default.minp*scale(j), img,...
                            [ Default.OutPutFolder,Default.filename,'.ppm'],Default.PpmOption) + 1;
        else
	   a = segmentImgOpt( Default.sigm*scale(j), Default.k*scale(j), Default.minp*scale(j), img,...
                            [ Default.OutPutFolder,Default.filename,'.ppm'], 0) + 1;
        end
        a = a(3:(end-2),3:(end-2)); %*** clean the edge superpixel index errors ***
        
        % Arrange the superpixel index in order
        if j == 1 % For the smallest Scale           

           ma = max(a(:));
           Unique_a = unique(a);
           SparseIndex = sparse(ma,1);
           SparseIndex(Unique_a) = 1:size(Unique_a);
           MedSup = full(SparseIndex(a));

           %Downsample to size as prediected depth map
           Sup{j} = imresize(MedSup,[Default.VertYNuDepth Default.HoriXNuDepth],'nearest');
           % clean superpixel section ====================================================================
           % merage all small and disconneted points in 1st scale segmentation
           [Sup{j},SupNeighborTableFlag] = premergAllsuperpixel_efficient(Sup{j}, Default);
           % ==============================================================
           % ===============================

        else  % o/w don't need the MedSup

           %Downsample to size size as prediected depth map
           a = imresize(a,[Default.VertYNuDepth Default.HoriXNuDepth],'nearest');
           ma = max(a(:));
           Unique_a = unique(a);
           SparseIndex = sparse(ma,1);
           SparseIndex(Unique_a) = 1:size(Unique_a);
           Sup{j} = full(SparseIndex(a));
        end
        clear a SparseIndex Unique_a ma;


        % show superpixel
        if DisplayFlag
           figure(1);
           imagesc(Sup{j});
           newmap = rand(max(max(Sup{j})),3);
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
        
end % end of for j=1:3   

return;
