% compare Gaussian smoothing to diffusion smoothing
img = imread('house.jpg'); img = rgb2gray(img);

% diffusion smoothing
num_iter = 2*[2 20 80 130]; delta_t = 1/7; kappa = 30; option = 1;
figure;
for i = 1:max(num_iter);
    img = anisodiff2D(img, 1, delta_t, kappa, option);
    [member, position] = ismember(i, num_iter);
    if member
        subplot(2, length(num_iter), position);
        imshow(img, []); box on;
    end
    i
end

% Gaussian smoothing
sigma = sqrt(2*num_iter); hsize = [55, 55];
img = imread('house.jpg'); img = rgb2gray(img); A = img;
for i = 1:length(sigma);
    h = fspecial('gaussian', hsize, sigma(i));
    A = imfilter(img, h);
    subplot(2, length(num_iter), i+length(num_iter));
    imshow(A, []); box on;
    i
end
