%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment 4: Stereo and Segmentation
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frame1L = imread('Materials/frame_1L.png');
frame1LR = imread('Materials/frame_1LR.png');
frame1R = imread('Materials/frame_1R.png');
frame1RL = imread('Materials/frame_1RL.png');

% Task 1. Calculate disparity using the SSD algorithm
window_size = 1;

% Matlab's disparity function
disparity_ml = disparity(frame1L, frame1R)

