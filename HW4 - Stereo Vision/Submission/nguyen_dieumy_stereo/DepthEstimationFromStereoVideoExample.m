%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all; clear all;

%% Depth Estimation From Stereo Video
% This example shows how to detect people in video taken with a calibrated stereo 
% camera and determine their distances from the camera.
%% Load the Parameters of the Stereo Camera
% Load the |stereoParameters| object, which is the result of calibrating the 
% camera using either the |stereoCameraCalibrator| app or the |estimateCameraParameters| 
% function. 

% Load the stereoParameters object.
load('handshakeStereoParams.mat');

% Visualize camera extrinsics.
% showExtrinsics(stereoParams);

%% Create Video File Readers and the Video Player
% Create System Objects for reading and displaying the video.

videoFileLeft = 'handshake_left.avi';
videoFileRight = 'handshake_right.avi';

readerLeft = vision.VideoFileReader(videoFileLeft, 'VideoOutputDataType', 'uint8');
readerRight = vision.VideoFileReader(videoFileRight, 'VideoOutputDataType', 'uint8');
% player = vision.DeployableVideoPlayer('Location', [20, 400]);

%% Read and Rectify Video Frames
% The frames from the left and the right cameras must be rectified in order 
% to compute disparity and reconstruct the 3-D scene. Rectified images have horizontal 
% epipolar lines, and are row-aligned. This simplifies the computation of disparity 
% by reducing the search space for matching points to one dimension. Rectified 
% images can also be combined into an anaglyph, which can be viewed using the 
% stereo red-cyan glasses to see the 3-D effect.

frameLeft = readerLeft.step();
frameRight = readerRight.step();

[frameLeftRect, frameRightRect] = ...
        rectifyStereoImages(frameLeft, frameRight, stereoParams);

% figure;
% imshow(stereoAnaglyph(frameLeftRect, frameRightRect));
% title('Rectified Video Frames');

%% Compute Disparity
% In rectified stereo images any pair of corresponding points are located on 
% the same pixel row. For each pixel in the left image compute the distance to 
% the corresponding pixel in the right image. This distance is called the disparity, 
% and it is proportional to the distance of the corresponding world point from 
% the camera.

frameLeftGray  = rgb2gray(frameLeftRect);
frameRightGray = rgb2gray(frameRightRect);

% Global variable(s) for all tests
max_disparity = 64;  % Disparity range: 0 - 64

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1. Calculate disparity using the SSD algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_sizes_ssd = [1, 5, 11];

% Create 2x2 subplot of various SSD disparity maps
figure;

% Matlab's built in disparity
subplot(2,2,1);
disparityMap = disparity(frameLeftGray, frameRightGray);
imshow(disparityMap, [0, 64]);
title('Disparity Map MATLAB');
colormap jet
colorbar
hold on;

% SSD disparity for varing window sizes
for plot_i = 2:4
    window_size = window_sizes_ssd(plot_i-1);
    subplot(2,2,plot_i);
    disparity_map_ssd = calculate_disparity_ssd(frameLeftGray, frameRightGray, ...
                                      window_size, max_disparity);
    imshow(disparity_map_ssd, [0, max_disparity]);
    formatSpec = 'Disparity Map SSD (%dx%d)';
    str = sprintf(formatSpec, window_size, window_size);
    title(str);
    colormap jet
    colorbar
end

% Save figure to png
saveas(gcf,'Results/SSD_disparity_maps.png')
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 2. Calculate disparity using the NCC algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_sizes_ncc = [3, 5, 7];

% Create 2x2 subplot of various SSD disparity maps
figure;

% Matlab's built in disparity
subplot(2,2,1);
disparityMap = disparity(frameLeftGray, frameRightGray);
imshow(disparityMap, [0, 64]);
title('Disparity Map MATLAB');
colormap jet
colorbar
hold on;

% NCC disparity for varing window sizes
for plot_i = 2:4
    window_size = window_sizes_ncc(plot_i-1);
    subplot(2,2,plot_i);
    disparity_map_ncc = calculate_disparity_ncc(frameLeftGray, frameRightGray, ...
                                      window_size, max_disparity);
    imshow(disparity_map_ncc, [0, max_disparity]);
    formatSpec = 'Disparity Map NCC (%dx%d)';
    str = sprintf(formatSpec, window_size, window_size);
    title(str);
    colormap jet
    colorbar
end

% Save figure to png
saveas(gcf,'Results/NCC_disparity_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 3. Uniqueness constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try: For a stereo match to be considered unique, 
% the minimum matching cost times a uniqueness factor must 
% be smaller than the cost for the next best match
% Source: https://nerian.com/products/stereo-vision-core/downloads/svc_data_sheet.pdf 

% Arbitrary values, probably in range 1-inf
uniqueness_factors = [1, 1.2, 1.5, 2];  
window_size = 5;

figure;
for plot_i = 1:length(uniqueness_factors)
    factor = uniqueness_factors(plot_i);
    disparity_map_uniqueness = calculate_disparity_uniqueness(frameLeftGray, frameRightGray, ...
                                               max_disparity, window_size, factor);
    subplot(2,2,plot_i);
    imshow(disparity_map_uniqueness, [0, max_disparity]);
    formatSpec = 'Disparity Map - Uniqueness factor: (%d)';
    str = sprintf(formatSpec, factor);
    title(str);
    colormap jet
    colorbar
end

% Save figure to png
saveas(gcf,'Results/Uniqueness_disparity_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 4. Disparity smoothness constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjacent points on line should have similar disparities
% Disparity values should change slowly 
% Try: Choose the disparity with lowest smoothness cost + data cost
% Smoothness cost = (if d_p = d_q = 0) (if d_p =/= d_q = 1)
% d and p are pairs of pixels in window/neighborhood
% Smoothness_term(d1, d2) indicates the cost of assigning 
% disparities d1 and d2 to neighboring pixels
% Smoothness threshold defines allowed max difference?
% Choose disparity map with lowest cost
% Source: https://www.cs.auckland.ac.nz/~rklette/CCV-CIMAT/pdfs/B20-StereoMatchingPart3.pdf
% Doesn't work too well

smoothness_thresholds = [5, 10, 20, 50];
window_size = 1;

figure;
for plot_i = 1:length(smoothness_thresholds)
    thres = smoothness_thresholds(plot_i);
    disparity_map_smoothness = calculate_disparity_smoothness(frameLeftGray, frameRightGray, ...
                                               max_disparity, window_size, thres);
    subplot(2,2,plot_i);
    imshow(disparity_map_smoothness, [0, max_disparity]);
    formatSpec = 'Disparity Map - Smoothness constraint: (%d)';
    str = sprintf(formatSpec, thres);
    title(str);
    colormap jet
    colorbar
end

% Save figure to png
saveas(gcf,'Results/Smoothness_disparity_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 5. Generate outliers map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine disparity by performing left-right 
% consistency check
window_size = 5;
disparity_1 = calculate_disparity_ssd(frameLeftGray, frameRightGray,...
                              window_size, max_disparity);
disparity_1_flipped = fliplr(calculate_disparity_ssd(fliplr(frameRightGray), fliplr(frameLeftGray), ...
                                           window_size, max_disparity));
disparity_2 = calculate_disparity_ncc(frameLeftGray, frameRightGray,...
                              window_size, max_disparity);
disparity_2_flipped = fliplr(calculate_disparity_ncc(fliplr(frameRightGray), fliplr(frameLeftGray), ...
                                                     window_size, max_disparity));
outlier_map1 = calculate_outlier(disparity_1, disparity_1_flipped, 1);
outlier_map2 = calculate_outlier(disparity_2, disparity_2_flipped, 1);
figure;
subplot(2,1,1)
imshow(outlier_map1);
title('Outlier map for SSD (5x5)');
subplot(2,1,2)
imshow(outlier_map2);
title('Outlier map for NCC (5x5)');

% Save figure to png
saveas(gcf,'Results/Outlier_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 6. Compute depth from disparity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
points3D_SSD = compute_depth(disparity_1, stereoParams);
points3D_NCC = compute_depth(disparity_2, stereoParams); 

figure;
subplot(2,1,1)
imshow(points3D_SSD);
title('Depth Matrix for SSD (5x5)');
subplot(2,1,2)
imshow(points3D_NCC);
title('Depth Matrix for NCC (5x5)');

% Save figure to png
saveas(gcf,'Results/Depth_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 7. Synthetic stereo sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test on synthetic stereo sequences and compare 
% with ground truth disparity map

% Read in synthetic left/right images
Ileft = (imread("Materials/Ileft.png"));
Iright = (imread("Materials/Iright.png"));

% Ground truth disparity map
ILR = imread("Materials/frame_1LR.png");

window_size = 15;

% 1. Generate disparity ssd and ncc maps, same window size
% Also use matlab built-in function to compare

% Matlab's disparity
figure;
subplot(3,1,1);
disparity_map_synthetic_ml = disparity(Ileft, Iright);
imshow(disparity_map_synthetic_ml, [0, max_disparity]);
title('Disparity Map MATLAB');
colormap jet
colorbar
hold on;

subplot(3,1,2);
disparity_map_synthetic_ssd = calculate_disparity_ssd(Ileft, Iright, window_size, disp_range);
imshow(disparity_map_synthetic_ssd, [0, max_disparity]);
title('Disparity Map SSD (15x15)');
colormap jet
colorbar

subplot(3,1,3);
disparity_map_synthetic_ncc = calculate_disparity_ncc(Ileft, Iright, window_size, disp_range);
imshow(disparity_map_synthetic_ncc, [0, max_disparity]);
title('Disparity Map NCC (15x15)');
colormap jet
colorbar
hold off;

% Save figure to png
saveas(gcf,'Results/Synthetic_disparity_maps.png')   

% 2. Calculate map of errors, compare against ground-truth
ssd_map = uint8(disparity_map_synthetic_ssd);
ncc_map = uint8(disparity_map_synthetic_ncc);

subplot(3,1,1);
imshow(ILR, []);
title('Ground truth');
colorbar
hold on;

subplot(3,1,2);
error_ssd = abs(ILR - ssd_map);
imshow(error_ssd, []);
title('SSD error map');
colorbar

subplot(3,1,3);
error_ncc = abs(ILR - ncc_map);
imshow(error_ncc, []);
title('NCC error map');
colorbar

% Save figure to png
saveas(gcf,'Results/Error_maps.png')   

% 3. Histogram of disparity differences
figure;
subplot(2,1,1);
ssd_diff = ILR - ssd_map;
histogram(ssd_diff);
title('SSD histogram of disparity differences');

subplot(2,1,2);
ncc_diff = ILR - ncc_map;
histogram(ncc_diff);
title('NCC histogram of disparity differences');

% Save figure to png
saveas(gcf,'Results/Disparity_histograms.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Programming 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 8A. Disparity matching along each epipolar line
% Task 8B. Backtracking
% Task 8C. Displaying the disparity map, again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window_size = 5;
occ = 0.01;

disparity_map = calculate_disparity_dp(Ileft, Iright, max_disparity, window_size, occ);
disparity_map = display_dmap(disparity_map);
figure;
imshow(disparity_map, [0, max_disparity]);
title('Disparity Map - Dynamic Programming');
colormap jet
colorbar

% Save figure to png
saveas(gcf,'Results/DP_disparity_maps.png')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Other Matlab code ----------------------%
% % Reconstruct the 3-D Scene
% Reconstruct the 3-D world coordinates of points corresponding to each pixel 
% from the disparity map.
% 
% points3D = reconstructScene(disparityMap, stereoParams);
% 
% Convert to meters and create a pointCloud object
% points3D = points3D ./ 1000;
% ptCloud = pointCloud(points3D, 'Color', frameLeftRect);
% 
% Create a streaming point cloud viewer
% player3D = pcplayer([-3, 3], [-3, 3], [0, 8], 'VerticalAxis', 'y', ...
%     'VerticalAxisDir', 'down');
% 
% Visualize the point cloud
% view(player3D, ptCloud);
% % Detect People in the Left Image
% Use the |vision.PeopleDetector| system object to detect people.
% 
% Create the people detector object. Limit the minimum object size for
% speed.
% peopleDetector = vision.PeopleDetector('MinSize', [166 83]);
% 
% Detect people.
% bboxes = peopleDetector.step(frameLeftGray);
% % Determine The Distance of Each Person to the Camera
% Find the 3-D world coordinates of the centroid of each detected person and 
% compute the distance from the centroid to the camera in meters.
% 
% Find the centroids of detected people.
% centroids = [round(bboxes(:, 1) + bboxes(:, 3) / 2), ...
%     round(bboxes(:, 2) + bboxes(:, 4) / 2)];
% 
% Find the 3-D world coordinates of the centroids.
% centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
% X = points3D(:, :, 1);
% Y = points3D(:, :, 2);
% Z = points3D(:, :, 3);
% centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];
% 
% Find the distances from the camera in meters.
% dists = sqrt(sum(centroids3D .^ 2));
%     
% Display the detected people and their distances.
% labels = cell(1, numel(dists));
% for i = 1:numel(dists)
%     labels{i} = sprintf('%0.2f meters', dists(i));
% end
% figure;
% imshow(insertObjectAnnotation(frameLeftRect, 'rectangle', bboxes, labels));
% title('Detected People');
% % Process the Rest of the Video
% Apply the steps described above to detect people and measure their distances 
% to the camera in every frame of the video.
% 
% while ~isDone(readerLeft) && ~isDone(readerRight)
%     Read the frames.
%     frameLeft = readerLeft.step();
%     frameRight = readerRight.step();
%     
%     Rectify the frames.
%     [frameLeftRect, frameRightRect] = ...
%         rectifyStereoImages(frameLeft, frameRight, stereoParams);
%     
%     Convert to grayscale.
%     frameLeftGray  = rgb2gray(frameLeftRect);
%     frameRightGray = rgb2gray(frameRightRect);
%     
%     Compute disparity. 
%     disparityMap = disparity(frameLeftGray, frameRightGray);
%     
%     Reconstruct 3-D scene.
%     points3D = reconstructScene(disparityMap, stereoParams);
%     points3D = points3D ./ 1000;
%     ptCloud = pointCloud(points3D, 'Color', frameLeftRect);
%     view(player3D, ptCloud);
%     
%     Detect people.
%     bboxes = peopleDetector.step(frameLeftGray);
%     
%     if ~isempty(bboxes)
%         Find the centroids of detected people.
%         centroids = [round(bboxes(:, 1) + bboxes(:, 3) / 2), ...
%             round(bboxes(:, 2) + bboxes(:, 4) / 2)];
%         
%         Find the 3-D world coordinates of the centroids.
%         centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
%         X = points3D(:, :, 1);
%         Y = points3D(:, :, 2);
%         Z = points3D(:, :, 3);
%         centroids3D = [X(centroidsIdx), Y(centroidsIdx), Z(centroidsIdx)];
%         
%         Find the distances from the camera in meters.
%         dists = sqrt(sum(centroids3D .^ 2, 2));
%         
%         Display the detect people and their distances.
%         labels = cell(1, numel(dists));
%         for i = 1:numel(dists)
%             labels{i} = sprintf('%0.2f meters', dists(i));
%         end
%         dispFrame = insertObjectAnnotation(frameLeftRect, 'rectangle', bboxes,...
%             labels);
%     else
%         dispFrame = frameLeftRect;
%     end
%     
%     Display the frame.
%     step(player, dispFrame);
% end
% 
% Clean up.
% reset(readerLeft);
% reset(readerRight);
% release(player);
% % Summary
% This example showed how to localize pedestrians in 3-D using a calibrated 
% stereo camera.
% % References
% [1] G. Bradski and A. Kaehler, "Learning OpenCV : Computer Vision with the 
% OpenCV Library," O'Reilly, Sebastopol, CA, 2008.
% 
% [2] Dalal, N. and Triggs, B., Histograms of Oriented Gradients for Human 
% Detection. CVPR 2005.
% 
% _Copyright 2013-2014 The MathWorks, Inc._