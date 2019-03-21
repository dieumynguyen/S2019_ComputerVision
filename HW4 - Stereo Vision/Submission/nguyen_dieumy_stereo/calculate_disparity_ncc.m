%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map = calculate_disparity_ncc(left_img, right_img, window_size, max_disparity)
        
    % Get size of left img as size of map & initialize disparity map
    rows = size(left_img, 1);
    cols = size(left_img, 2);
    disparity_map = zeros(rows, cols);
    
    % Convert to double
    left_img = im2double(left_img);
    right_img = im2double(right_img);
        
    % Pad img
    padding = NaN(rows, max_disparity);
    right_img = horzcat(right_img, padding);
    
    % Center of window
    pix_center = (window_size - 1) / 2;
    
    % Loop over img rows and cols
    for row = 1+pix_center:rows-pix_center
        for col = 1+pix_center:cols-pix_center 
            % Left 
            left = left_img(row-pix_center:row+pix_center, col-pix_center:col+pix_center);
            left_mean = mean2(left);
            left_norm = left - left_mean;
            left_sd = sqrt(sum(sum((left - left_mean).^2)));
            if left_sd == 0, left_sd = 1; end
            left = repmat(left_norm / left_sd, 1, max_disparity+1);
            
            % Right
            right = right_img(row-pix_center:row+pix_center, col-pix_center:col+pix_center);
            right_mean = mean2(right);
            right_norm = right - right_mean;
            right_sd = sqrt(sum(sum((right - left_mean).^2)));
            if right_sd == 0, right_sd = 1; end
            right = right_norm / right_sd;

            % Calculate NCC
            for d = 1:max_disparity
               c = col + d;
               window = right_img(row-pix_center:row+pix_center, c-pix_center:c+pix_center);
               window_mean = mean2(window);
               window_norm = window - window_mean;
               window_sd = sqrt(sum(sum((window - window_mean).^2)));
               if window_sd == 0, window_sd = 1; end
               window = window_norm / window_sd;
               right = horzcat(right, window);
            end
            NCCs = sum(reshape(left.*right, window_size.^2, max_disparity+1));
            
            % Find best/highest disparity
            disparity = find(NCCs == max(NCCs));
            disparity_map(row, col) = disparity(1) - 1;
        end
    end  
end
