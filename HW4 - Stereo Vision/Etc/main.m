%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment 4: Stereo and Segmentation
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

%----------- Read images
I1 = imread('Materials/frame_1L.png');
I2 = imread('Materials/frame_1R.png');

%----------- Convert to grayscale
I1gray = rgb2gray(I1);
I2gray = rgb2gray(I2);

%----------- Rectify images
blobs1 = detectSURFFeatures(I1gray, 'MetricThreshold', 2000);
blobs2 = detectSURFFeatures(I2gray, 'MetricThreshold', 2000);

[features1, validBlobs1] = extractFeatures(I1gray, blobs1);
[features2, validBlobs2] = extractFeatures(I2gray, blobs2);

indexPairs = matchFeatures(features1, features2, 'Metric', 'SAD', ...
                            'MatchThreshold', 5);

matchedPoints1 = validBlobs1(indexPairs(:,1),:);
matchedPoints2 = validBlobs2(indexPairs(:,2),:);

[fMatrix, epipolarInliers, status] = estimateFundamentalMatrix(...
  matchedPoints1, matchedPoints2, 'Method', 'RANSAC', ...
  'NumTrials', 10000, 'DistanceThreshold', 0.1, 'Confidence', 99.99);

if status ~= 0 || isEpipoleInImage(fMatrix, size(I1)) ...
  || isEpipoleInImage(fMatrix', size(I2))
  error(['Either not enough matching points were found or '...
         'the epipoles are inside the images. You may need to '...
         'inspect and improve the quality of detected features ',...
         'and/or improve the quality of your images.']);
end

inlierPoints1 = matchedPoints1(epipolarInliers, :);
inlierPoints2 = matchedPoints2(epipolarInliers, :);

[t1, t2] = estimateUncalibratedRectification(fMatrix, ...
  inlierPoints1.Location, inlierPoints2.Location, size(I2));
tform1 = projective2d(t1);
tform2 = projective2d(t2);

[I1Rect, I2Rect] = rectifyStereoImages(I1, I2, tform1, tform2);

% ----------- Task 1. Calculate disparity using the SSD algorithm
I1RectGray = rgb2gray(I1Rect);
I2RectGray = rgb2gray(I2Rect);

window_sizes = [1, 5, 11];

disparity_map = SDD_disparity(I1RectGray, I2RectGray, 100, 0);
figure; imshow(disparity_map);
colormap(gca,jet) 
colorbar

% % Matlab's disparity function
% disparity_ml = disparity(I1RectGray, I2RectGray);
% figure; imshow(disparity_ml); title('Disparity Map MATLAB');
% colormap(gca,jet) 
% colorbar
