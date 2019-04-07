%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = ComputeFeatures(img)
    % Compute a feature vector for all pixels of an image. You can use this
    % function as a starting point to implement your own custom feature
    % vectors.
    %
    % INPUT
    % img - Array of image data of size h x w x 3.
    %
    % OUTPUT
    % features - Array of computed features for all pixels of size h x w x f
    %            such that features(i, j, :) is the feature vector (of
    %            dimension f) for the pixel img(i, j, :).

    img = double(img);
    height = size(img, 1);
    width = size(img, 2);
    features = zeros([height, width, 8]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    %                              YOUR CODE HERE                             %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    features(:, :, 1:5) = ComputePositionColorFeatures(img);
    features(:, :, 6:7) = ComputeGradientFeatures(img);
    features(:, :, 8) = ComputeEdgeFeatures(img); 
end
