%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 3 Question 8
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write a script that compares the performance of a 3x3 mean
% filter, a 3x3 median filter and Gaussian filters with different 
% values of ?. Get familiar with the following MATLAB functions: 
% imnoise, medfilt2, conv2, filter2, fspecial, imfilter, and edge.

% clear all; close all

% Load image and convert to grayscale
img = imread("peppers.png");
img_grayscale = rgb2gray(img);

% Add salt & pepper noise to image
% subplot(3,2,1);
% img_noise = imnoise(img_grayscale, 'salt & pepper');
% imshow(img_noise);
% title("Salt & pepper image");

% Add Gaussian white noise with mean 0 and ? = 1 
% in the [0, 255] range (or ? = 1/256 in the [0, 1] range) 
img_double = double(img_grayscale) / 255;
subplot(3,2,1);
img_noise = imnoise(img_double, 'gaussian', 0, 1/256);
imshow(img_noise);
title("Gaussian noise image");

% Mean filter:
subplot(3,2,2);
img_mean = filter2(fspecial('average', 3), img_noise) / 255;
imshow(img_mean);
title("Mean filter");

% Median filter:
subplot(3,2,3);
img_median = medfilt2(img_noise);
imshow(img_median);
title("Median filter");

% Gaussian filters:
% sigma = 1/3 pixel
subplot(3,2,4);
H_1 = fspecial('gaussian', 3, 1/3);
img_gauss_1 = imfilter(img_noise, H_1,'replicate');
imshow(img_gauss_1);
title("Gaussian filter, \sigma=1/3");

% sigma = 1 pixel
subplot(3,2,5);
H_2 = fspecial('gaussian', 3, 1);
img_gauss_2 = imfilter(img_noise, H_2,'replicate');
imshow(img_gauss_2);
title("Gaussian filter, \sigma=1");

% sigma = 1.5 pixel
subplot(3,2,6);
H_3 = fspecial('gaussian', 3, 1.5);
img_gauss_3 = imfilter(img_noise, H_3,'replicate');
imshow(img_gauss_3);
title("Gaussian filter, \sigma=1.5");

