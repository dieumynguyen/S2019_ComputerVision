%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = ComputeColorFeatures(img)
    % Computes a feature vector of color values for all pixels in an image.
    %
    % INPUTS
    % img - Array of image data of size h x w x 3.
    %
    % OUTPUT
    % features - Array of size h x w x 3 where features(i, j, :) is the feature
    %            vector for the pixel img(i, j, :). In the case of color
    %            features this is simply the input image.

    features = double(img);
end