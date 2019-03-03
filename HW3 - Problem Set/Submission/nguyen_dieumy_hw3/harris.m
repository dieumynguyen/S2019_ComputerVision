%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 3 Question 9
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Harris corner detection algorithm
function [corner_coords, descriptors] = harris(I, w, threshold, suppression)
    % 1. Filter with Gaussian window, W, for initial de-noising
    H = fspecial('gaussian', w);
    bw = imfilter(I, H);

    % 2. Compute magnitude of gradient everywhere
    dy = [-1 0 1; -1 0 1; -1 0 1];
    dx = dy';
                 
    % Compute Ix and Iy image derivatives by convolving
    Ix = conv2(bw, dx, 'same');
    Iy = conv2(bw, dy, 'same');

    % 3. Compute M matrix using the same W Gaussian window function
    Ix2 = conv2(Ix .^ 2, H, 'same'); 
    Iy2 = conv2(Iy .^ 2, H, 'same');
    Ixy = conv2(Ix .* Iy, H, 'same');

    % 4. Compute the cornerness
    cornerness = (Ix2 .* Iy2 - Ixy.^2) ./ (Ix2 + Iy2 + eps);

    % 5. Find points whose surrounding window gave large corner response
    % & 6. Take points of local maxima (non-max suppression)
    if suppression == 0    % If no non-max suppression, only threshold
        cornerness = cornerness > threshold;
        [rows, columns] = find(cornerness);
    elseif suppression == 1 % Else apply non-max suppression with threshold
        order = 3;
        im_ordfilt = ordfilt2(cornerness, order^2, ones(order)); % Dilate image
        cornerness = (cornerness == im_ordfilt) & (cornerness > threshold);
        [rows, columns] = find(cornerness);
    end
    
    % Store rows and columns in output matrix 
    corner_coords = [rows columns];
    
    % Make matrix for pixel intensities of 3x3 patch around each point
    [imsize_x, imsize_y] = size(I);
    for i = 1:length(corner_coords)
        x = corner_coords(i,1);
        y = corner_coords(i,2);
        window_x = x-1:x+1;
        window_x(window_x <= 0) = window_x(2);
        window_x(window_x > imsize_x) = window_x(2);
        window_y = y-1:y+1;
        window_y(window_y <= 0) = window_y(2);
        window_y(window_y > imsize_y) = window_y(2);
        descriptor = reshape(bw(window_x, window_y), 1, []);
        descriptors(i, :) = descriptor * 256; % Back to grayscale range
    end
end
