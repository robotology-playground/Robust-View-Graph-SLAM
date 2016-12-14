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
function [relativeFeatureVector] = makeRelativeFeatureVector(H,scales)

global GeneralDataFolder ScratchDataFolder LocalFolder ClusterExecutionDirectory...
    ImgFolder VertYNuPatch VertYNuDepth HoriXNuPatch HoriXNuDepth a_default b_default Ox_default Oy_default...
    Horizon_default filename batchSize NuRow_default SegVertYSize SegHoriXSize WeiBatchSize PopUpVertY PopUpHoriX taskName;

%global columnWidth rowWidth 
global nDim nLaw
global nStatistics
global nHistBins
global minHist maxHist stepHist

if nargin < 2
	scale = 1;
end
%if nargin < 2
%	type = 'laws';
%end

nLaw = size(H,3);
numscales = length(scales);
H1size = nLaw;%nDim / nStatistics;

nHistBins = 10;
%load([GeneralDataFolder '/maxHist.mat']);
%load([GeneralDataFolder '/minHist.mat']);
load([GeneralDataFolder '/MM.mat']);
if scales == 1
	maxHist = Max{1}';
	minHist = Min{1}';
elseif scales == 2
	maxHist = Max{2}';
	minHist = Min{2}';
else
	maxHist = Max{3}';
	minHist = Min{3}';
end
maxHist = log(1+maxHist);
minHist = log(1+minHist);
maxHist = log(1+maxHist);
minHist = log(1+minHist);

%maxHist = [17, 12, 12, 12, 12, 12, 12, 12, 12,    15.9, 15.9,       24, 24, 24, 24, 24, 24]'; % How to pick values
%minHist = [10,  0,  0,  0,  0,  0,  0,  0,  0,     14.7, 14.7,        8,  8,  8,  8,  8,  8]';% How to pick values
maxHist = maxHist.^2;
minHist = minHist.^2;
stepHist = (maxHist - minHist) / (nHistBins-1);
step = 1/nHistBins;
      

%HistMinMax;
relativeFeatureVector = zeros( VertYNuDepth, HoriXNuDepth, H1size*numscales*nHistBins );

%for s=1:numscales
        s = 1;
	scale = scales(1);
	%reductionScale = 1/(2*(scale-1) + 1);
	%if scale ~= 1
	%	resizedImg = imresize(img, reductionScale, 'bilinear');
	%else
%		resizedImg = img;
%	end
	%resizedImg(:,:,1) = medfilt2(resizedImg(:,:,1), [5, 5], 'symmetric');
	%resizedImg(:,:,2) = medfilt2(resizedImg(:,:,2), [5, 5], 'symmetric');
	%resizedImg(:,:,3) = medfilt2(resizedImg(:,:,3), [5, 5], 'symmetric');

%	numoverlaps = 0;
%	if scale == 2
%		numoverlaps = scale - 1;
%	end
	numoverlaps = scale - 1;

	edgefactor = (2*numoverlaps + 1);

	% Assume that the image is correctly oriented.
	imheight = size(H, 1);
	imwidth  = size(H, 2);
	% Step sizes in x and y
	stepwidth  = imwidth/HoriXNuDepth;
	stepheight = imheight/VertYNuDepth;

	H = log( 1  +  H  );
	H = log( 1  +  H  );
   for l=1:nLaw
    H(:,:,l) = (H(:,:,l) - repmat( minHist(l), size(H,1), size(H,2)) ) ./ ...
                        repmat( maxHist(l)-minHist(l), size(H,1), size(H,2) );
   end
	%====================================================================
	%==================The Laws' Filters Applied=========================
	%====================================================================

	intstepheight = floor((2*numoverlaps+1)*stepheight);
	intstepwidth  = floor((2*numoverlaps+1)*stepwidth);
	for g = 1:VertYNuDepth
    	for c = 1:HoriXNuDepth
			gridLeft  = round( (g - 1 - numoverlaps)*stepheight + 1);
			gridRight = round( (g + numoverlaps)*stepheight );
			tempstep = gridRight - gridLeft + 1;
			residue = tempstep - intstepheight;
			if residue == 1
				gridRight = gridRight - residue;
			end
			if residue > 1 || residue < 0
				residue
				display('oops');
			end
			normfactor = 1.0;
			if gridLeft < 1
				gridLeft = 1;
				normfactor = normfactor * edgefactor / (numoverlaps + g);
			elseif gridRight > (VertYNuDepth * stepheight)
				gridRight = round(VertYNuDepth * stepheight);
				normfactor = normfactor * edgefactor ...
					/ (numoverlaps + VertYNuDepth - g + 1);
			end
			gridTop = round( (c - 1 - numoverlaps)*stepwidth + 1);
			gridBot = round( (c + numoverlaps)*stepwidth );
			tempstep = gridBot - gridTop + 1;
			residue = tempstep - intstepwidth;
			if residue == 1
				gridBot = gridBot - residue;
			end
			if residue > 1 || residue < 0
				residue
				display('oops');
			end
			if gridTop < 1
				gridTop = 1;
				normfactor = normfactor * edgefactor / (numoverlaps + c);
			elseif gridBot > (HoriXNuDepth * stepwidth)
				gridBot = round(HoriXNuDepth * stepwidth);
				normfactor = normfactor * edgefactor ...
					/ (numoverlaps + HoriXNuDepth - c + 1);
            end

            
      		% The Laws Histogram
			tmpPatch = reshape( abs(H(gridLeft:gridRight, gridTop:gridBot, :)), ...
                                [(gridRight-gridLeft+1)*(gridBot-gridTop+1) nLaw]);
                   %         histc( tmpPatch, min
            rfv = histc( tmpPatch, [-inf, step:step:(1-step), inf]);
            rfv = permute( reshape( rfv(1:nHistBins,:), nHistBins*H1size, 1), [2 3 1]);
    		relativeFeatureVector(g,c,...
					((s-1)*nHistBins*H1size + 1):( (s-1)*nHistBins*H1size + H1size*nHistBins) ) = ...
					rfv;
		end
	end
	clear H;

%end
return;
