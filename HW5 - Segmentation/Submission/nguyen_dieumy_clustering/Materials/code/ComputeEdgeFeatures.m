%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = ComputeEdgeDetectionFeatures(img)
    height = size(img, 1);
    width = size(img, 2);
    
    features = zeros([height, width, 1]);
    gray_img = rgb2gray(img);
 
    features(:,:,1) = edge(gray_img, 'Canny');
end