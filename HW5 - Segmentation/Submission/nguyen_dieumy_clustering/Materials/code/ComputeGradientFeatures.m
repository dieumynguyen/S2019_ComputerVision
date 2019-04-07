%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = ComputeGradientFeatures(img)
    height = size(img, 1);
    width = size(img, 2);
    
    gray_img = rgb2gray(img);
    features = zeros([height, width, 2]);
    
    [mag, dir] = imgradient(gray_img);
    features(:,:,1) = mag;
    features(:,:,2) = dir;
end