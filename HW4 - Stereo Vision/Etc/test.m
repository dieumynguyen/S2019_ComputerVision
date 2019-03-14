% Perform Basic Block Matching.
% clear all
% close all

% Read images.
leftI = frameLeftGray;
rightI = frameRightGray;
figure(1), imshow(leftI, []), title('Left image');
figure(2), imshow(rightI, []), title('Right image');

% Convert to floating point.
leftI = single(leftI);
rightI = single(rightI);


% For every pixel in the right image, we extract the 7-by-7-pixel block
% around it and search along the same row in the left image for the block
% that best matches it. Here we search in a range of +-disparityRange
% pixels around the pixel's location in the first image, and we use the sum
% of absolute differences (SAD) to compare the image regions. We need only
% search over columns and not over rows because the images are rectified.
% We use the |TemplateMatcher| System object to perform this block matching
% between each block and the region of interest.

D = zeros(size(leftI), 'single');      % Disparity map

disparityRange = 64;

% Blocksize is (2*halfBlockSize+1)-by-(2*halfBlockSize+1).
halfBlockSize = 3;
blockSize = 2*halfBlockSize+1;

% Allocate space for all template matcher System objects.
tmats = cell(blockSize);

% Initialize progress bar
hWaitBar = waitbar(0, 'Performing basic block matching...');
nRowsLeft = size(leftI, 1);

% Scan over all rows.
for m=1:nRowsLeft
    % Set min/max row bounds for image block.
    minr = max(1,m-halfBlockSize);
    maxr = min(nRowsLeft,m+halfBlockSize);
    % Scan over all columns.
    for n=1:size(leftI,2)
        minc = max(1,n-halfBlockSize);
        maxc = min(size(leftI,2),n+halfBlockSize);
        % Compute disparity bounds.
        mind = max( -disparityRange, 1-minc );
        maxd = min( disparityRange, size(leftI,2)-maxc );

        % Construct template and region of interest.
        template = rightI(minr:maxr,minc:maxc);
        templateCenter = floor((size(template)+1)/2);
        roi = [minc+templateCenter(2)+mind-1 ...
               minr+templateCenter(1)-1 ...
               maxd-mind+1 1];
        % Lookup proper TemplateMatcher object; create if empty.
        if isempty(tmats{size(template,1),size(template,2)})
            tmats{size(template,1),size(template,2)} = ...
                vision.TemplateMatcher('ROIInputPort',true);
        end
        thisTemplateMatcher = tmats{size(template,1),size(template,2)};
        
        % Run TemplateMatcher object.
        loc = step(thisTemplateMatcher, leftI, template, roi);
        D(m,n) = loc(1) - roi(1) + mind;
    end
    waitbar(m/nRowsLeft,hWaitBar);
end
close(hWaitBar);


% Show resulting disparity map in pseudocolor.
% For display purposes, we saturate the depth map to have only positive
% values. In general, slight angular misalignment of the stereo cameras
% used for image acquisition can allow both positive and negative
% disparities to appear validly in the depth map. In this case, however,
% the stereo cameras were near perfectly parallel, so the true disparities
% have only one sign. Thus this correction is valid.
figure
imshow(D,[]), axis image, colormap('jet'), colorbar, impixelinfo;
caxis([0 disparityRange]);
title('Depth map from basic block matching');