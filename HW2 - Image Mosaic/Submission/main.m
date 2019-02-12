%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(42);

clear all; close all; clc;

% Global variables
num_points = 10;
num_trials = 20;

% Images from different views
img1 = imread('Materials/Square0.jpg');
img2 = imread('Materials/Square1.jpg');

%------ Task 1. Getting correspondences
% Uncomment 2 lines below to select new points
% coordinates_matrix_1 = getPoint(img1, img2, num_points);
% save("coordinates_matrix_1.mat", "coordinates_matrix_1");

% Else use saved points
load('coordinates_matrix_1.mat');

%------ Task 2. Computing homography
H = computeH(coordinates_matrix_1, num_points, num_trials);

%------ Task 3. Warping between image planes
[warped_img, bb_xmin, bb_ymin] = warp1(img1, H);

%------ Task 4. Create the output mosaic
mosaic_img = mosaic(img2, warped_img, bb_xmin, bb_ymin, 0, 0);

%------ Task 5.1. Display original images and mosaic
figure;
subplot(2,2,1);
imagesc(img1);
title('Image 1');

subplot(2,2,2);
imagesc(img2);
title('Image 2');

subplot(2,2,[3,4]);
imagesc(mosaic_img);
title('Squared Mosaic');
imwrite(mosaic_img, 'Results/SquaredMosaic.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------ Task 5.2. Two more mosaics

%%%%% Image set 1: Bee %%%%%
img1_bee = imread('Materials/Bee1.jpg');
img2_bee = imread('Materials/Bee2.jpg');
% Uncomment 2 lines below to select new points
% coordinates_matrix_2 = getPoint(img1_bee, img2_bee, num_points);
% save("coordinates_matrix_2.mat", "coordinates_matrix_2");
load('coordinates_matrix_2.mat');
H_bee = computeH(coordinates_matrix_2, num_points, num_trials);
[warped_img_bee, bb_xmin_bee, bb_ymin_bee] = warp1(img1_bee, H_bee);
mosaic_img_bee = mosaic(img2_bee, warped_img_bee, bb_xmin_bee, bb_ymin_bee, 50, 25);
figure;
subplot(2,2,1);
imagesc(img1_bee);
title('Image 1');
subplot(2,2,2);
imagesc(img2_bee);
title('Image 2');
subplot(2,2,[3,4]);
imagesc(mosaic_img_bee);
title('Bee Mosaic');
imwrite(mosaic_img_bee, 'Results/BeeMosaic.jpg');

%%%%% Image set 2: Hike %%%%%
img1_hike = imread('Materials/Hike1.jpg');
img2_hike = imread('Materials/Hike2.jpg');
% Uncomment 2 lines below to select new points
% coordinates_matrix_3 = getPoint(img1_hike, img2_hike, num_points);
% save("coordinates_matrix_3.mat", "coordinates_matrix_3");
load('coordinates_matrix_3.mat');
H_hike = computeH(coordinates_matrix_3, num_points, num_trials);
% Using warp2() and mosaic2(0 functions which are not optimized
% but work better for these images than warp1() and mosaic()
[warped_img_hike, up, left] = warp2(img1_hike, H_hike);
[mosaic_img_hike] = mosaic2(img2_hike, warped_img_hike, up, left);
figure;
subplot(2,2,1);
imagesc(img1_hike);
title('Image 1');
subplot(2,2,2);
imagesc(img2_hike);
title('Image 2');
subplot(2,2,[3,4]);
imagesc(mosaic_img_hike);
title('Hike Mosaic');
imwrite(mosaic_img_hike, 'Results/HikeMosaic.jpg');

%------ Task 5.3. Warp 1 image into a frame of the other
% Photo credits: Bartush, Photography Pro
num_points = 4;
img1_frame = imread('Materials/Frame1.jpg');
img2_frame = imread('Materials/Frame2.jpg');
% Uncomment 2 lines below to select new points
% coordinates_matrix_4 = getPoint(img1_frame, img2_frame, num_points);
% save("coordinates_matrix_4.mat", "coordinates_matrix_4");
load('coordinates_matrix_4.mat');
H_frame = computeH(coordinates_matrix_4, num_points, num_trials);
[warped_img_frame, up, left] = warp2(img1_frame, H_frame);
[mosaic_img_frame] = mosaic2(img2_frame, warped_img_frame, up, left);
figure;
subplot(2,2,1);
imagesc(img1_frame);
title('Image 1');
subplot(2,2,2);
imagesc(img2_frame);
title('Image 2');
subplot(2,2,[3,4]);
imagesc(mosaic_img_frame);
title('Frame Mosaic');
imwrite(mosaic_img_frame, 'Results/FrameMosaic.jpg');

