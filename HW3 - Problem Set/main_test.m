
clear all; close all;

% Algorithm variables
w = 3;
sigma = 2;
thresh = 0.1;
sze = 11;
disp = 0;

% Load and convert image to black and white
im = imread('checkerboard.png');
bw = double(im(:,:,1)) ./ 256;

% 1. Filter with Gaussian window, W, for initial de-noising
H = fspecial('gaussian', w);
bw = imfilter(bw, H);

% 2. Compute magnitude of gradient everywhere
% Define the filters to compute the image derivatives
dy = [-1 0 1; -1 0 1; -1 0 1];
dx = dy'; % dx is the transpose matrix of dy

% Image derivatives Ix and Iy are the horizontal and vertical edges
% Compute Ix and Iy by convolving 
Ix = conv2(bw, dx, 'same');
Iy = conv2(bw, dy, 'same');

% 3. Compute M matrix for each image window with same W Gaussian
% window function
% Calculating the gradient of the image Ix and Iy
% g = fspecial('gaussian', max(1,fix(6*sigma)), sigma);
Ix2 = conv2(Ix .^ 2, H, 'same'); % Smoothed squared image derivatives
Iy2 = conv2(Iy .^ 2, H, 'same');
Ixy = conv2(Ix .* Iy, H, 'same'); 

% 4. Compute the cornerness measure (R) for each image window
cornerness = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);

% 5 & 6. Non-maximal suppression (take points of local max)
% and threshold (find where R > threshold)
mx = ordfilt2(cornerness, sze^2, ones(sze)); % Grey-scale dilate
cornerness = (cornerness == mx) & (cornerness > thresh); % Find maxima
[rws, cols] = find(cornerness); % Find row,col coords.

clf; imshow(bw);
hold on;
p=[cols rws];
plot(p(:,1),p(:,2),'gs',...
    'MarkerSize',8,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor', 'r');
title('\bf Harris Corners')