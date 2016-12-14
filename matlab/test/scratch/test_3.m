

% Notes:
% (1) xx1 and xx2 should both be of size n-by-2 matrices, with one image point per row.
%     Here n is the number of corresponding points and, for the estimation of the
%     fundamental matrix, n >= 7.  In the code for rectification, we follow the
%     convention that each column contains an image point in homogeneous coordinates.
% (2) xx1 and xx2 should be in x-y coordinates, with the origin of the image coordinate
%     system at the top-left corner and y-axis pointing down.  For the computation
%     that follows, we need to convert the image coordinate system to be origined
%     at the centre of the image buffer with the y-axis pointing down.











%image/epipolar lines rectification

% matches
inliers=false(1,length(matches));
inliers(inlier_matches)=true;
pts1=kps1(matches(1,inliers),1:2);
pts2=kps2(matches(2,inliers),1:2);
hpts1 = makehomogeneous(pts1');
hpts2 = makehomogeneous(pts2');







% %% rectification code 5: using vl-feat applications
% this projets the right image into the left image using homography
% [rect_im1,rect_im2,H]=rect_right_to_left(im1,im2,pts1,pts2);
% figure; imshow((rect_im1+rect_im2)/2);


% % left image
% figure;
% subplot(2,2,1);
% image(im1); axis image; title('left'); hold on;
% plot(pts1(:,1), pts1(:,2),'r+','MarkerSize',5);
% 
% % right image
% subplot(2,2,2);
% image(im2); axis image; title('right'); hold on;
% plot(pts2(:,1), pts2(:,2),'g+','MarkerSize',5);
% 
% % Epipolar geometry
% [e1,e2]=epipoles_from_F(FM); % (FM' * e1 ~ 0) and (FM  * e2 ~ 0)
% e1_=e1./e1(3);
% for i=1:size(pts1,1)
%     % left image
%     subplot(2,2,1);
%     left_x=pts1(i,1);left_y=pts1(i,2);
%     left_p=[left_x;left_y;1];
%     left_epipolar_x=1:2*size(im2,2);
%     % using the eqn of line: ax+by+c=0; y = (-c-ax)/b
%     left_epipolar_y=left_y+(left_epipolar_x-left_x)*...
%         (e1_(2)-left_y)/(e1_(1)-left_x);
%     plot(left_epipolar_x,left_epipolar_y,'r');
%     % right image
%     subplot(2,2,2);
%     right_p=FM*left_p;
%     right_epipolar_x=1:2*size(im2,2);
%     right_epipolar_y=(-right_p(3)-right_p(1)*right_epipolar_x)/right_p(2);
%     plot(right_epipolar_x,right_epipolar_y,'g');
% end
% 
% % rectifying transformation
% cy = round( rows/2 );
% cx = round( cols/2 );
% [H1,H2]=epipolar_rectification(FM,pts1,pts2,cx,cy);
% 
% % rectify the images
% Img1_new = rectify_image(im1,H1);
% Img2_new = rectify_image(im2,H2);
% 
% % plot
% figure;
% subplot(1,2,1);imagesc(Img1_new);axis image;grid on;
% subplot(1,2,2);imagesc(Img2_new);axis image;grid on;
% 
% % plot rectified points
% hpts1 = makehomogeneous(pts1');
% hpts1 = H1*hpts1;
% npts1 = hnormalise(hpts1);
% hpts2 = makehomogeneous(pts2');
% hpts2 = H2*hpts2;
% npts2 = hnormalise(hpts2);
% figure;
% plot(npts1(1,:),npts1(2,:),'r+');
% hold on; grid on;
% plot(npts2(1,:),npts2(2,:),'go');
% 
% return;
% 
% %% plot the epipolar lines of randomly selected points in the image
% [m,n]=size(im1g);
% figure(1),imagesc(im1g); colormap(gray);
% title('Click a point on this Left Image');axis image;
% figure(2),imagesc(im2g); colormap(gray);
% title('Corresponding Epipolar Line in this Right Image');axis image;
% list=['r' 'b' 'g' 'y' 'm' 'k' 'w' 'c'];
% figure(1); hold on;
% left_x=zeros(size(list,2),1);
% left_y=zeros(size(list,2),1);
% for i=1:size(list,2);
%     % clicking a point on the left image:
%     [left_x(i),left_y(i)]=ginput(1);
%     plot(left_x(i),left_y(i),'r*');
%     
%     % finding the epipolar line on the right image:
%     left_p=[left_x(i);left_y(i);1];
%     right_p=FM*left_p; 
%     right_p=right_p/right_p(3);
%     right_epipolar_x=1:2*m;
%  
%     % using the eqn of line: ax+by+c=0; y = (-c-ax)/b
%     right_epipolar_y=(-right_p(3)-right_p(1)*right_epipolar_x)/right_p(2);
%     figure(2);
%     hold on;
%     plot(right_epipolar_x,right_epipolar_y,list(mod(i,8)+1));
%     
%     % now finding the other epipolar line we know that left epipole is the
%     % 3rd column of V. we get V from svd of F. F=UDV'
%     left_epipole=FV(:,3);
%     left_epipole=left_epipole/left_epipole(3);
%     left_epipolar_x=1:2*m;
%     left_epipolar_y=left_y(i)+(left_epipolar_x-left_x(i))*...
%         (left_epipole(2)-left_y(i))/(left_epipole(1)-left_x(i));
%     figure(1);
%     plot(left_epipolar_x,left_epipolar_y,list(mod(i,8)+1));
% end
% 
% 
