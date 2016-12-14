% plot support points
%a=load('support_pt.txt');
%figure;
%subplot(2,2,1); imagesc(im1g); title('Left image');
%hold on; plot(a(:,1),a(:,3),'m+');
%subplot(2,2,2); imagesc(im2g); title('right image');
%hold on; plot(a(:,2),a(:,3),'m+');
%subplot(2,2,3); imagesc(D1'); title('Left disparity');
%subplot(2,2,4); imagesc(D2'); title('right disparity');

% plot estimated disparity
[X,Y]=meshgrid(1:size(im1,1),1:size(im1,2));
X=X(1:1:end)/20; Y=Y(1:1:end)/20; D=D1(1:1:end);
ind=D>0; X=X(ind); Y=Y(ind); D=D(ind);

% apply arbitrary rotation (due to unknown camera orientation)
r=0*pi/180; p=0*pi/180; q=0*pi/180;
R=angle2dcm(r,p,q,'XYZ');
points=R*[X; Y; D];
X=points(1,:);
Y=points(2,:);
Z=points(3,:);

% plot with no colors
figure;
plot3(D,Y,-X,'.','markersize',.5);
xlabel('Z-axis');ylabel('Y-axis');zlabel('-X-axis');
axis equal; grid on;

% plot with colors
depth=sqrt(X.^2+Y.^2+D.^2);
[sorted, I] = sort(depth);
sorted_X = X(I);
sorted_Y = Y(I);
sorted_D = D(I);
temp_min = min(min(sorted));
step = 1;
numColors=150;

color2 = hsv(numColors); % color by depth/height
figure; hold on
for iii = 1:numColors
        temp_index = find(sorted >= temp_min & sorted < (temp_min + step));
        temp_X = sorted_X(temp_index);
        temp_Y = sorted_Y(temp_index);
        temp_D = sorted_D(temp_index);
        temp_min = temp_min + step;
        plot3(temp_D,temp_Y,-temp_X,'.','markersize',.5,'color',color2(iii,:));
end
xlabel('Z-axis');ylabel('Y-axis');zlabel('-X-axis');
axis equal; grid on;
