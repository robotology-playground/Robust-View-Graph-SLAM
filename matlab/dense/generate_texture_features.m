function [TextureFeature, TextSup] = generate_texture_features(param, im, SupMedScale, SmallSup, HiSup, FeatureType)

% This function calculate all features that using the
% output from calculateFilterBanks_old (Time and ram comsuming)
% These Feature including AbsFeature DiffFeature HistogramFeature
% FeatureTypr:            1          2           3

% Interface:
% Input---
% 1) Default is the structure contain all the default parameter setting
% 2) Img is the RGB image matrix
% 3) SmallSup is Superpixel index matrix of size [VertYNuDepth,HoriXNuDepth]
% 4) HiSup is high resolution Superpixel index matrix of size [SegVertYSize,SegHoriXSize] or smaller
% 5) FeatureType AbsFeature(1) DiffFeature(2) HistogramFeature(3)
% Output---
% 1) TextureFeature.Abs:
% 2) TextureFeature.Diff:
% 3) TextureFeature.Hist:

% Flag set up for testing purpose ========================================
% Must be stablized after testing XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Default.Flag.NormalizeFlag = 1; % Default.Flags.NormalizeFlag set to 1
% if we want to normalize the Feature according to the Trainning Set.
DominateSupFilter = 1; % DominateSupFilter set to 1 if we want to use only
% the dominate superpixel area for acerage feature calculation.
% ========================================================================

% load FeaMax to do normalizeing
load([param.ParaFolder 'FeaMax.mat']);
if param.Flag.NormalizeFlag == 1
    FeaMax = 10.^floor(log10(FeaMax));
else
    FeaMax(:) = 1;
end

% Initial Setting of TextureFeature Structure:
% size(TextureFeature.Abs) = (No of Depth Point, No of Features: 17 texture filters, H2 and H4 kinds, 3 Scales)
TextureFeature.Abs = zeros(param.VertYNuDepth*param.HoriXNuDepth, 17*2*3);
TextureFeature.Diff = zeros(1,1);
TextureFeature.Hist = zeros(1,1);

% change the resolution in to exactly [TrainVerYSize TrainHoriXSize]
if ~all(size(im) == [param.TrainVerYSize,param.TrainHoriXSize,3])
    im = imresize(im,[param.TrainVerYSize,param.TrainHoriXSize],'bilinear');
end

% Start the loop for 3 Scale
fInd = 2;
for i = 1:3
    
    % only if DominateSupFilter and the largest scale will use superpixel mask
    if DominateSupFilter && i==1
        SupMaskFlag = 1;
    else
        SupMaskFlag = 0;
    end
    
    % Generate the 17 texture filter output
    % set to global for InnerMultiple sup generation
    global H2;
    H2 = calculateFilterBanks_old(im); % (hard work 1min) use Ashutaosh's code
    % InnerMulSup
    if i == 2
        TextSup = gen_TextSup_efficient_GlobalH2(param);
    end
    
    % Scale Setting
    ResY = round(param.TrainVerYSize/(3^(i))); % ResY ResX will be the new image size
    ResX = round(param.TrainHoriXSize/(3^(i)));
    im = imresize(im,[ResY,ResX],'nearest');
    
    % Now H2 stand for H2
    H2 = H2.^2;
    if any(FeatureType == 1) % calculate the AbsFeatures
        StartAbs = 34*(i-1)+1;
        [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)), fInd] = ...
            AbsFeatureGenMex_MemoryEfficient(param, SmallSup, HiSup, SupMaskFlag, FeaMax, fInd);
    end
    if any(FeatureType == 3) % calculate the HistogramFeature
        relativeFeatureVector = makeRelativeFeatureVector(H2, 1);
    end
    if any(FeatureType == 2) % calculate the DiffFeature
    end
    
    % This H2 is H4
    H2 = H2.^2;
    if any(FeatureType == 1) % calculate the AbsFeatures
        StartAbs = 34*(i-1)+1+17;
        [TextureFeature.Abs(:, (StartAbs):(StartAbs+16)), fInd] = AbsFeatureGenMex_MemoryEfficient(...
            param, SmallSup, HiSup, SupMaskFlag, FeaMax, fInd);
    end
    
    clear global H2;
    
end

% GoundBoundaryFea
RayCenter = GenerateRay(param.HoriXNuDepth, param.VertYNuDepth, 'center',...
    param.a_default, param.b_default, param.Ox_default, param.Oy_default); %[horiXSizeLowREs,VertYSizeLowREs,3]
big_sup = imresize(SupMedScale, [param.VertYNuDepth, param.HoriXNuDepth], 'nearest');
DefaultGroundMask = [zeros(floor(param.VertYNuDepth/2), param.HoriXNuDepth); ...
    ones(param.VertYNuDepth-floor(param.VertYNuDepth/2), param.HoriXNuDepth)];
GroundSupIndex = unique(big_sup(logical(DefaultGroundMask)));
% GroundSupIndex = unique(big_sup(logical(maskg)));
for j = 1:param.HoriXNuDepth
    Gmask = logical(zeros(size(big_sup(:,j))));
    for k = GroundSupIndex'
        Gmask(big_sup(:,j) == k) = true;
    end
    [rowSub, colSub ] = find(Gmask);
    minRow = min(rowSub);
    
    GroundVertEdge(:,j) = (1-(param.VertYNuDepth-1)/param.VertYNuDepth)./RayCenter(:,j,3);
    if size(minRow,1) ~= 0
        GroundVertEdge((minRow+1):param.VertYNuDepth,j) = (1-((minRow+1):param.VertYNuDepth)'...
            /param.VertYNuDepth)./RayCenter((minRow+1):param.VertYNuDepth,j,3);
        GroundVertEdge(1:(minRow),j) = (1- minRow/param.VertYNuDepth)./RayCenter(1:(minRow),j,3);
    end
end
GroundVertEdgePics = GroundVertEdge;

TextureFeature.Abs = [TextureFeature.Abs min([GroundVertEdgePics(:)...
    repmat((param.VertYNuDepth:-1:1)'/param.VertYNuDepth,[param.HoriXNuDepth 1])./...
    reshape(RayCenter(:,:,3),[],1)],[],2)];

TextureFeature.Abs = [SmallSup(:) TextureFeature.Abs];
