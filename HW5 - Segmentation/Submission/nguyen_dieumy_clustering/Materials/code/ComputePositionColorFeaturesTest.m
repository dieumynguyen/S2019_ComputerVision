%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ComputePositionColorFeaturesTest()
    % Tests the implementation of ComputePositionColorFeatures by comparing the
    % current implementation with the output from the reference implementation.

    load('../test_data/ComputePositionColorFeaturesTest.mat');
    
    your_features = ComputePositionColorFeatures(img);
    if all(your_features == features)
        disp(['Congrats! Your implementation of ComputePositionColorFeatures ' ...
              'appears to work correctly.']);
    else
        disp(['Uh oh - Your implementation of ComputePositionColorFeatures ' ...
              'does not produce the same output as the reference implementation.']);
    end

end