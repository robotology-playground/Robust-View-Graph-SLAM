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
function []=gen_Sup_new(sigm,k,min,SelectSegmentationPara);
% this function generate superpixel using default parameter
% but can also change to manually input parameter

%%% Jeff's Comments
% The global variable filename specifies a vector of names.  For
% each file specified, read in the jpg and write out a ppm of size
% SegVertYSize x SegHoriXSize.  Send that image to the CMU
% segmentation program with params 0.8*[sigm, k, min].  If
% SelectSegmenationPara is true, then display the results and ask
% the user to enter new sigm, k, and min; repeat until user is happy.
%
% Resize both the image and the CMU output to PopUpVertY x
% PopUpHoriX and write them to the scratch/ppm folder.
%
% Take the unresized CMU output and call suprgb2ind, converting it
% from an image to a matrix of superpixel indicies.  save this 
% index image to scratch/data/MedSeg/ and a low res (VertYNuDepth x
% HoriXNuDepth) version to scratch/data/.
%%%%


% default parameter
if nargin < 4
    SelectSegmentationPara = 0; % if SelectSegmentationPara == 1, enable the parameter interation with user.
end

% declaim global variable
global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize PopUpVertY PopUpHoriX;

scale =[0.8 1.6 5]; % use different scale to generate small(0.8) middle(1.6) 5(large) scale of superpixel

% generate superpixel of each image

NuPics = size(filename,2);
for i = 1:NuPics
%for i = 1:10
%    i  
    for j = 1:1%3% number of scale of superpixel

        
%         sigm_new =
        % load image and process it to Hi Medi and Low Resolution
        Img = imread([GeneralDataFolder '/' ImgFolder '/' filename{i} '.jpg']); % Readin the high resolution image
        [VertYSizeHiREs HoriXSizeHiREs dummy]= size(Img);% find the dimension size of the Hi Resolution image
        clear dummy;
        % Loadin the GroundTruth data to know the depthMap size
%         depthfile = strrep(filename{i},'img','depth'); % the depth filename(without .file extension) associate with the *jpg file
%         load([GeneralDataFolder '/depthMap/' depthfile '.mat']); 
%         [VertYSizeLowREs HoriXSizeLowREs]= size(depthMap);% find the dimension size of the depth data
        % in the new laser data we have scatter depthmap so use a
        % predecided LowRes
        VertYSizeLowREs = VertYNuDepth;
        HoriXSizeLowREs = HoriXNuDepth;
        
        % using a fixed range of median size image [SegVertYSize SegHoriXSize ] 
        %  to generate superpixel to reduce computation
        if VertYSizeHiREs*HoriXSizeHiREs > SegVertYSize*SegHoriXSize
%             Ratio = [VertYSizeHiREs HoriXSizeHiREs]./[SegVertYSize SegHoriXSize]
%              Ratio = [1 1]
%             [VertYSizeHiREs HoriXSizeHiREs ]./floor(Ratio)
%             pause
            Img = imresize(Img,[SegVertYSize+1 SegHoriXSize+1 ],'nearest'); % Downsample high resolution image to a fixed median size image
            %Img = imresize(Img,[VertYSizeHiREs HoriXSizeHiREs ]./floor(Ratio),'nearest'); % Downsample high resolution image to a range of median size image
            imwrite(Img,[ScratchDataFolder '/ppm/' filename{i} '.ppm'],'ppm');% store median Resolution image to PPM format to feed in CMU C++ function
        else
            imwrite(Img,[ScratchDataFolder '/ppm/' filename{i} '.ppm'],'ppm');% store median Resolution image to PPM format to feed in CMU C++ function
        end    
        
        % choose superpixel of the images
        % default segmentation parameter
        ok = 0; % ok ==1 means accept the segmentation
        while 1
            % call segment function writen in C++ from MIT
            disp([LocalFolder '/../third_party/Superpixels/segment ' num2str(sigm*scale(j)) ' ' num2str(k*scale(j)) ...
                ' ' num2str(min*scale(j)) ' ' ScratchDataFolder '/ppm/' filename{i} '.ppm' ' ' ...
                ScratchDataFolder '/ppm/' filename{i} '_' num2str(sigm*scale(j)) '_' ...
                num2str(k*scale(j)) '_' num2str(min*scale(j)) '.ppm']);
            system([LocalFolder '/../third_party/Superpixels/segment ' num2str(sigm*scale(j)) ' ' num2str(k*scale(j)) ...
                ' ' num2str(min*scale(j)) ' ' ScratchDataFolder '/ppm/' filename{i} '.ppm' ' ' ...
                ScratchDataFolder '/ppm/' filename{i} '_' num2str(sigm*scale(j)) '_' ...
                num2str(k*scale(j)) '_' num2str(min*scale(j)) '.ppm']);
%            system([LocalFolder '/../third_party/Superpixels/segment ' num2str(sigm*scale(j)) ' ' num2str(k*scale(j)) ...
%                ' ' num2str(min*scale(j)) ' ' ScratchDataFolder '/ppm/' filename{i} '.ppm' ' ' ...
%                ScratchDataFolder '/ppm/' filename{i} '.ppm']);

            MediResImgSuperpixel = imread([ScratchDataFolder '/ppm/' filename{i} '_' num2str(sigm*scale(j)) '_' num2str(k*scale(j)) '_' num2str(min*scale(j)) '.ppm']); % Readin the high resolution image
            MediResImgSuperpixel = MediResImgSuperpixel(1:(end-1),1:(end-1),1:3);
            figure(1); image(MediResImgSuperpixel); % show the superpixel in Medi Resolution
                imwrite(MediResImgSuperpixel,[ScratchDataFolder '/ppm/' 'test' int2str(i) '.jpg'],'jpg');
    
            % check if need to select segmentation parameter
            if SelectSegmentationPara==1;
                ok = input('Is the segmentation of image OK');% input new segmentation parameter
            else    
                ok =1 ;% accept default segmentation parameter
            end
    
            % finish segmentation clean up the ppm folder.
            if ok==1;
                delete([ScratchDataFolder '/ppm/' filename{i} '_' num2str(sigm*scale(j)) '_' num2str(k*scale(j)) '_' num2str(min*scale(j)) '.ppm']);
%                delete([ScratchDataFolder '/ppm/' filename{i} '.ppm']);
                
                % generating [PopUpVertY, PopUpHoriX] (default 800 x 600) 
                % .ppm superpixel and image for photopopup software.
                newfilename = strrep(filename{i},'.','')
                SupPopup = imresize(MediResImgSuperpixel,[PopUpVertY, PopUpHoriX],'nearest');
                imwrite(SupPopup,[ScratchDataFolder '/ppm/' newfilename '.ppm'],'ppm');
                ImgPopup = imresize(Img,[PopUpVertY, PopUpHoriX],'nearest');
                imwrite(ImgPopup,[ScratchDataFolder '/ppm/' newfilename '.jpg'],'jpg');
%            return;    
                break;
            end
            sigm = input('type sigm of segmentation');
            k = input('type k of segmentation');
            min = input('type min of segmentation');
            
        end
% ============used to test no residual superpixel================
%    end;
%end;    
%return;
% ===============================================================

        % index superpixel
        [MediResImgIndexSuperpixelSepTemp dummy]= suprgb2ind(MediResImgSuperpixel); clear dummy;
        LowResImgIndexSuperpixelSep{i,j} = imresize(MediResImgIndexSuperpixelSepTemp,...
            [VertYSizeLowREs HoriXSizeLowREs],'nearest'); %Downsample to size size as prediected depth map
        
        % merage all small point in higher scale segmentation
        if j ~= 1
           LowResImgIndexSuperpixelSep{i,j} = premergAllsuperpixel(LowResImgIndexSuperpixelSep{i,j});
        end
        if j == 1;
	   MediResImgIndexSuperpixelSep = MediResImgIndexSuperpixelSepTemp;
        end
        % refining superpixel
        % superpixel segmentation LowResImgSeperatedSuperpixel
        %LowResImgsuperpixel = imresize(MediResImgSuperpixel,[VertYSizeLowREs HoriXSizeLowREs],'nearest'); %Downsample high resolution image to the same pixel size of GroundTruth data
        %[LowResImgIndexSuperpixel  LowResImgIndexSuperpixel_list]= suprgb2ind(LowResImgsuperpixel);

        % comment: cmu's superpixel might be connected. use premergsuperpixel to
        % deal with nonconnected superpixels and very small superpixels
        %[LowResImgIndexSuperpixelSepTemp]=premergsuperpixel(LowResImgIndexSuperpixel); % hard work 1min
        
        % reorder the index number of the LowResImgIndexSuperpixelSep
        %[LowResImgIndexSuperpixelSep{i,j}  LowResImgIndexSuperpixelSep_list]= ordersup(LowResImgIndexSuperpixelSepTemp);
        
        % show superpixel
        figure(2);
        imagesc(LowResImgIndexSuperpixelSep{i,j});
        newmap = rand(max(max(LowResImgIndexSuperpixelSep{i,j})),3);
        colormap(newmap);
        
        % process the MediResImgSuperpixel to have the same number of
        % LowResImgIndexSuperpixelSep
%         if j==1 
%             tic
%             [MediResImgIndexSuperpixel dummy]= suprgb2ind(MediResImgSuperpixel); clear dummy;
%             MediResImgIndexSuperpixelSep = imresize(LowResImgIndexSuperpixelSep{i,1},size(MediResImgIndexSuperpixel),'nearest');
%             NuSupMedi = max(max(MediResImgIndexSuperpixel));
%             LowToMediResImgIndexSuperpixel = zeros(size(MediResImgIndexSuperpixel));
%             for k = 1:NuSupMedi
%                 mask = MediResImgIndexSuperpixel==k;
%                 LowToMediResImgIndexSuperpixel(mask) = analysesupinpatch(MediResImgIndexSuperpixelSep(mask));
% %                 [list_sup] = analysesupinpatch(MediResImgIndexSuperpixelSep(mask));
% %                 [I C] = max(list_sup(2,:));
% %                 LowToMediResImgIndexSuperpixel(mask) = list_sup(1,C);
%             end
%             LowToMediResImgIndexSuperpixelSep{i} =...
%             premergAllsuperpixel(LowToMediResImgIndexSuperpixel);
%             toc
%         end    
% ===============
    end    
save([ScratchDataFolder '/data/MedSeg/MediResImgIndexSuperpixelSep' num2str(i) '.mat'], 'MediResImgIndexSuperpixelSep');
end    
% =============

% save result for later application
save([ScratchDataFolder '/data/LowResImgIndexSuperpixelSep.mat'], 'LowResImgIndexSuperpixelSep');

return;
