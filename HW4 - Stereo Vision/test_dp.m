close all, clear all

% Read in images
left_I = imread('Materials/Ileft.png');
right_I = imread('Materials/Iright.png');
left_I = mean(left_I, 3);
right_I = mean(right_I, 3);

% Algorithm variables
disp_range = 64;
h_block_size = 3; 
penalty = 0.01;

I_disp = stereo_dp(left_I, right_I, disp_range, h_block_size, penalty);
I_disp(I_disp < 70) = nan;
figure; imshow(display_dmap(I_disp), [0 disp_range/3]); colormap('jet'); colorbar;
title('Disparity map via dynamic programming')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

