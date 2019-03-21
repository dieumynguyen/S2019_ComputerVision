%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map = calculate_disparity_uniqueness(left_img, right_img, max_disparity, ... 
                                                        window_size, uniqueness_factor)
                                                    
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
           left = left_img(row-pix_center:row+pix_center, col-pix_center:col+pix_center);
           left = repmat(left, 1, max_disparity+1);
           right = right_img(row-pix_center:row + pix_center, col-pix_center:col+pix_center);
           for d = 1:max_disparity
              c = col+d;
              window = right_img(row+pix_center:row+pix_center, c-pix_center:c+pix_center);
              right = horzcat(right, window);
           end
           ssd = (left - right).^2;
           ssd = imgaussfilt(ssd, 'FilterSize', window_size);
           ssd = reshape(ssd, window_size.^2, max_disparity+1);

           % Sum disparities & find min
           ssd = sum(ssd);
           disparity = find(ssd == min(ssd));

           % Next best match
           sorted_ssd = sort(ssd(:));
           disparity_2 = sorted_ssd(2);

           % Uniqueness constraint: For a match to be unique,
           % the matching cost (SSD) times a uniqueness factor 
           % must be smaller than the cost for the next best match
           % This is like how much the best/min cost value should
           % "win" over the second best to be considered the correct
           % match
           if min(ssd) * uniqueness_factor < disparity_2
                best = disp(1);
                disparity_map(row, col) = best - 1;
           else
                disparity_map = disparity_map;
           end
       end
    end
end
