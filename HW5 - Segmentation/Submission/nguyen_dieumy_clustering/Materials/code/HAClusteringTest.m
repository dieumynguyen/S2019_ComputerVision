%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 5: Segmentation via Clustering
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HAClusteringTest(visualize)
% Tests your implementation of HAClustering.m by comparing its output on a
% test dataset with the output of our reference implementation on the same
% dataset.
%
% INPUT
% visualize - Whether or not to visualize each step of the clustering
%             algorithm. Optional; default is true.

    if nargin < 1
        visualize = true;
    end
    load('../test_data/HAClusteringTest.mat');
    my_idx = HAClustering(X, k, visualize);
    if all(my_idx == idx)
        disp(['Congrats! Your HAClustering algorithm produces the same ' ...
             'output as ours.']);
    else
        disp(['Uh oh - Your HAClustering algorithm produces a different ' ...
              ' output from ours.']);
    end
end