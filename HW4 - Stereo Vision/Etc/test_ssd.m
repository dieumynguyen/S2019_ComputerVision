% left=imread('Ileft.png');
% right=imread('Iright.png');

frameLeftGray  = rgb2gray(frameLeftRect);
frameRightGray = rgb2gray(frameRightRect);

left = imresize(frameLeftGray, 0.5);
right = imresize(frameRightGray, 0.5);

% Algorithm variables
disp_range = 64;
h_block_size = 3;

I_disp = stereo_ssd(left, right, h_block_size, disp_range);
imshow(I_disp, [0 64]);
colormap jet;
colorbar ;
title('Depth map from  block matching');
