%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map = calculate_disparity_ssd(left_img, right_img, window_size, max_disparity)
    
    % Initialize disparity map & get size of left img as size of map
    disparity_map = zeros(size(left_img), 'single');

    rows = size(left_img, 1);
    cols = size(left_img, 2);

    left_img = mean(left_img, 3);
    right_img = mean(right_img, 3);
    
    % Loop over img rows and cols
    for row = 1:rows
        row_max = max(1, row - window_size);
        row_min = min(rows, row + window_size);

        for col = 1:cols
              col_max = max(1, col - window_size);
              col_min = min(cols, col + window_size);

              % Pixel search limit
              pix_min = max(-max_disparity, 1 - col_max);
              pix_max = min(max_disparity, cols - col_min);

              template = right_img(row_max:row_min, col_max:col_min);

              window_count = pix_max - pix_min + 1;
              SSDs = zeros(window_count, 1);

              % Calculate squared difference
              for i = pix_min:pix_max
                  window = left_img(row_max:row_min, (col_max+i):(col_min+i));
                  idx = i - pix_min + 1;
                  SSDs(idx,1)= sumsqr(template - window);
              end
              
              % Find best disparity
              [~, ssd] = sort(SSDs);
              match_index = ssd(1,1);
              disparity = match_index + pix_min - 1;
              disparity_map(row, col) = disparity;
        end
    end
end