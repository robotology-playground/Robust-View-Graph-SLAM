function [kp1,kp2]=test_fast(im01,im02,options)

% make image greyscale
if length(size(im01))==3
    %im1=im01(:,:,2);
    im1=rgb2gray(im01);
else
    im1=im2double(im01);
end
if length(size(im02))==3
    %im2=im02(:,:,2);
    im2=rgb2gray(im02);
else
    im2=im2double(im02);
end

kp1 = []; kp2 = [];

% fast features
kp1=fast9(im1,options.fastthreshold,options.fastnonmax);
kp2=fast9(im2,options.fastthreshold,options.fastnonmax);
%POINTS = detectFASTFeatures(im1,'MinQuality',.01,'MinContrast',.001)
%POINTS = detectFASTFeatures(im2,'MinQuality',.01,'MinContrast',.001)
kp1=kp1(~isnan(kp1(:,1))&~isnan(kp1(:,2)),1:2);
kp2=kp2(~isnan(kp2(:,1))&~isnan(kp2(:,2)),1:2);

% image grid
if options.gridsize>0;
    xlimits=[options.gridmargin size(im1,2)-options.gridmargin];
    ylimits=[options.gridhorizon size(im1,1)-options.gridmargin];
    x=xlimits(1):options.gridsize:xlimits(2);
    y=ylimits(1):options.gridsize:ylimits(2);
    n=length(x)*length(y)+size(kp1,1);
    while n>options.maxnumfeat;
        options.gridsize=1.05*options.gridsize;
        x=xlimits(1):options.gridsize:xlimits(2);
        y=ylimits(1):options.gridsize:ylimits(2);
        n=length(x)*length(y)+size(kp1,1);
    end
    fprintf(['grid size changed to : ', num2str(options.gridsize),'\n']);
    [x,y]=meshgrid(x,y);
    kp1=[kp1;x(:),y(:)];
    kp2=[kp2;x(:),y(:)];
end

% % texture analysis using Local range of image
% K1 = rangefilt(im1);
% [y, x] = find(K1>20);
% kp1 = [kp1; x(1:3:end), y(1:3:end)];
% K2 = rangefilt(im2);
% [y, x] = find(K2>20);
% kp2 = [kp2; x(1:3:end), y(1:3:end)];

% % texture analysys using Local entropy of grayscale image
% J1 = entropyfilt(im1);
% [y, x] = find(J1>4.5);
% kp1 = [kp1; x(1:3:end), y(1:3:end)];
% J2 = entropyfilt(im2);
% [y, x] = find(J2>4.5);
% kp2 = [kp2; x(1:3:end), y(1:3:end)];

% % superpixels analysis
% superpixels_analysis(im01);

% limit scan to image margins
xlimits=[options.fastmargin size(im1,2)-options.fastmargin];
ylimits=[options.gridhorizon size(im1,1)-options.fastmargin];
kp1=kp1(kp1(:,1)>xlimits(1)&kp1(:,1)<xlimits(2)&...
   kp1(:,2)>ylimits(1)&kp1(:,2)<ylimits(2),:);
kp2=kp2(kp2(:,1)>xlimits(1)&kp2(:,1)<xlimits(2)&...
   kp2(:,2)>ylimits(1)&kp2(:,2)<ylimits(2),:);
fprintf([' ',num2str(size(kp1,1)),', ']);
fprintf([' ',num2str(size(kp2,1)),'\n']);

% show features and images
if options.verbose>1
    clf;
    subplot(1,2,1);image(im01);hold on;
    plot(kp1(:,1),kp1(:,2),'+');
    axis([1 size(im1,2) 1 size(im1,1)]);axis('ij');
    subplot(1,2,2);image(im02);hold on;
    plot(kp2(:,1),kp2(:,2),'+');
    axis([1 size(im2,2) 1 size(im2,1)]);axis('ij');
    drawnow;
end