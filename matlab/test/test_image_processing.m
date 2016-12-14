

% image processing
function [im1_, im2_, im1g_, im2g_] = test_image_processing(im1, im2)

% image channles
R1 = im1(:,:,1); R1_hist = imhist(R1);
G1 = im1(:,:,2); G1_hist = imhist(G1);
B1 = im1(:,:,3); B1_hist = imhist(B1);
R2 = im2(:,:,1); R2_hist = imhist(R2);
G2 = im2(:,:,2); G2_hist = imhist(G2);
B2 = im2(:,:,3); B2_hist = imhist(B2);

% figure;
% subplot(1, 3, 1); bar(R1_hist, 'r');
% subplot(1, 3, 2); bar(G1_hist, 'g');
% subplot(1, 3, 3); bar(B1_hist, 'b');
% figure;
% subplot(1, 3, 1); bar(R2_hist, 'r');
% subplot(1, 3, 2); bar(G2_hist, 'g');
% subplot(1, 3, 3); bar(B2_hist, 'b');

% hist equalisation
R1_ = histeq(R1); R1_hist_ = imhist(R1_);
G1_ = histeq(G1); G1_hist_ = imhist(G1_);
B1_ = histeq(B1); B1_hist_ = imhist(B1_);
R2_ = histeq(R2); R2_hist_ = imhist(R2_);
G2_ = histeq(G2); G2_hist_ = imhist(G2_);
B2_ = histeq(B2); B2_hist_ = imhist(B2_);

% figure;
% subplot(1, 3, 1); bar(R1_hist_, 'r');
% subplot(1, 3, 2); bar(G1_hist_, 'g');
% subplot(1, 3, 3); bar(B1_hist_, 'b');
% figure;
% subplot(1, 3, 1); bar(R2_hist_, 'r');
% subplot(1, 3, 2); bar(G2_hist_, 'g');
% subplot(1, 3, 3); bar(B2_hist_, 'b');

% new image
im1_(:,:,1) = R1_;
im1_(:,:,2) = G1_;
im1_(:,:,3) = B1_;
im2_(:,:,1) = R2_;
im2_(:,:,2) = G2_;
im2_(:,:,3) = B2_;

% convert images
im1g_ = rgb2gray(im1);
im2g_ = rgb2gray(im2);