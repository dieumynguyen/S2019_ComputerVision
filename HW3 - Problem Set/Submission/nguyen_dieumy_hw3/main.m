%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 3 Question 9 
% Algorithm: Harris corner detector
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

% Algorithm variables
w = 3;
threshold = 0.1; 
suppression = 1;

%--------- Part A. Testing on original image
im_name = "checkerboard.png";
im_orig = imread(im_name) + 10;
figure;
subplot(1,3,1); imshow(im_orig); title("Original image");
subplot(1,3,2); plot_harris(im_orig, w, threshold, 0);
subplot(1,3,3); plot_harris(im_orig, w, threshold, 1);

%--------- Part B. Testing on rotated, translated, scaled image
% 1. Rotation
im_rotated = imrotate(im_orig, 45);
figure;
subplot(1,2,1); plot_harris(im_rotated, w, threshold, 0);
subplot(1,2,2); plot_harris(im_rotated, w, threshold, 1);

% 2. Translated
im_translate = imtranslate(im_orig, [20 40]);
figure;
subplot(1,2,1); plot_harris(im_translate, w, threshold, 0);
subplot(1,2,2); plot_harris(im_translate, w, threshold, 1);

% 3. Scaled 
im_scaled = imresize(im_orig, 0.5);
figure;
subplot(1,2,1); plot_harris(im_scaled, w, threshold, 0);
subplot(1,2,2); plot_harris(im_scaled, w, threshold, 1);

%--------- Part C. Testing on 4 new versions
% 1. Brighter: Add constant positive offset to all pixels
im_bright_add = im_orig + 100;
figure;
subplot(1,2,1); plot_harris(im_bright_add, w, threshold, 0);
subplot(1,2,2); plot_harris(im_bright_add, w, threshold, 1);

% 2. Darker: Add constant negative offset to all pixels
im_dark_add = im_orig + -100;
figure;
subplot(1,2,1); plot_harris(im_dark_add, w, threshold, 0);
subplot(1,2,2); plot_harris(im_dark_add, w, threshold, 1);

% 3. Brighter: Multiply constant positive offset to all pixels
im_bright_mult = im_orig * 10;
figure;
subplot(1,2,1); plot_harris(im_bright_mult, w, threshold, 0);
subplot(1,2,2); plot_harris(im_bright_mult, w, threshold, 1);

% 4. Darker: Multiply constant negative offset to all pixels
im_dark_mult = im_orig * 0.5;
figure;
subplot(1,2,1); plot_harris(im_dark_mult, w, threshold, 0);
subplot(1,2,2); plot_harris(im_dark_mult, w, threshold, 1);






