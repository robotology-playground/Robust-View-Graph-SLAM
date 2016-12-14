function [Sup,MedSup,FeatureSup,TextureFeature] = superpixels_analysis(img, K, options)

% parameters
param.SegVertYSize = options.imgsize(1);
param.SegHoriXSize = options.imgsize(2);
param.TrainVerYSize = options.imgsize(1);
param.TrainHoriXSize = options.imgsize(2);
param.VertYNuDepth = options.imgsize(1)/2;%55;
param.HoriXNuDepth = options.imgsize(2)/2;%305;
param.VertYNuPatch = options.imgsize(1)/10;%55;
param.HoriXNuPatch = options.imgsize(2)/10;%61;
param.sigm = 0.3;%(.5,.8,.3)
param.k = 100;%(100,200,300)
param.minp = 100;%(100,150,20)
param.OutPutFolder = [pwd,'/'];
param.filename = 'superpixels';
param.PpmOption = 0;
param.SmallThre = 5;
param.scale = [.5,1,5];
param.ParaFolder = '/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/dense/';
param.Flag.NormalizeFlag = 0;
param.display = 0;
param.Flag.DisplayFlag = 0;
param.fx = K(1,1);
param.fy = K(2,2);
param.Ox = K(1,3);
param.Oy = K(2,3);
param.a_default = options.imgsize(1)/param.fy;%0.70783777; %0.129; % horizontal physical size of image plane normalized to focal length (in meter)
param.b_default = options.imgsize(2)/param.fx;%0.946584169;%0.085; % vertical physical size of image plane normalized to focal length (in meter)
param.Ox_default = 1-param.Ox/options.imgsize(1);%0.489272914; % camera origin offset from the image center in horizontal direction
param.Oy_default = 1-param.Oy/options.imgsize(2);%0.488886982; % camera origin offset from the image center in vertical direction
param.GroundThreshold = 0.5;
param.SkyThreshold = 1;

% 1) Basic Superpixel generation and Sup clean
%fprintf('Creating Superpixels...           ');
%[MedSup,Sup,param,SupNeighborTable] = gen_super_pixels(param,img);
[MedSup,Sup,param,SupNeighborTable] = gen_Sup_efficient(param, img);
%disp([ num2str( toc(startTime) ) ' seconds.']);

% 2) Texture Features and inner multiple Sups generation
%fprintf('Creating Features and multiple segmentations... ');
% [TextureFeature,TextSup] = generate_texture_features(param,img,Sup{2},Sup{1},...
%     imresize((MedSup),[param.TrainVerYSize,param.TrainHoriXSize],'nearest'), 1);
[TextureFeature,TextSup] = GenTextureFeature_InnerMulSup(param, img, Sup{2}, Sup{1},...
    imresize((MedSup),[param.TrainVerYSize param.TrainHoriXSize],'nearest'), 1);%, maskg);
%disp([ num2str( toc(startTime) ) ' seconds.']);

% 3) Superpixel Features generation
%fprintf('Calculating superpixel-shape features...       ');
FeatureSup = f_sup_old(param, Sup{1}, MedSup, SupNeighborTable);
%disp([ num2str( toc(startTime) ) ' seconds.']);

%if Default.Flag.IntermediateStorage
%   save([ ScratchFolder '/' strrep( filename{1},'.jpg','') '_IM.mat' ],'FeatureSup','TextureFeature','Sup','TextSup');
%end

%clf;
%subplot(2,2,1); imshow(Sup{1}./max(max(Sup{1})));
%subplot(2,2,2); imshow(Sup{2}./max(max(Sup{2})));
%subplot(2,2,3); imshow(Sup{3}./max(max(Sup{3})));
%subplot(2,2,4); imshow(MedSup./max(max(MedSup)));
%drawnow;