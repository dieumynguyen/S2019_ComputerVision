%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the disparity map between 2 images using DP with backtracking
function disparity_map = calculate_disparity_dp(left_img, right_img, max_disparity, window_size, penalty)
    
    % Initialize disparity map & get size of img
    disparity_map = zeros(size(left_img), 'single');

    rows = size(left_img, 1);
    cols = size(left_img, 2);

    left_img = mean(left_img, 3);
    right_img = mean(right_img, 3);
    
    % Keep track of cost
    diff_min = Inf;
    cost = diff_min * ones(cols, 2*max_disparity+1, 'single');

    % Compute disparity using dynamic programming with backtracking
    % Loop over rows
    for row = 1:rows
        cost(:) = diff_min;
        row_max = max(1, row-window_size);
        row_min = min(rows, row+window_size);

        % Loop over columns
        for col = 1:cols
            col_max = max(1, col-window_size);
            col_min = min(cols, col+window_size);

            % Set the pixel search limit with disparity range
            pix_min = max(-max_disparity, 1-col_max);
            pix_max = min(max_disparity, cols-col_min);

            % Search within that range, calc costs
            for px = pix_min:pix_max
                left_pix = left_img(row_max:row_min, (col_max:col_min)+px);
                right_pix = right_img(row_max:row_min, col_max:col_min);
                
                % Squared error measure bw pixels
                cost(col, px+max_disparity+1) = sumsqr(left_pix-right_pix);
            end
        end

        % Compute minimal cost alignment of 2 scanlines with backtracking 
        map = zeros(size(cost), 'single');
        cost_p = cost(end, :);
        for col = cols-1:-1:1
            diff = (cols - col + 1) * diff_min;
            [z, idx] = min([diff diff cost_p(1:end-4) + 3 * penalty;
                          diff cost_p(1:end-3) + 2 * penalty;
                          cost_p(1:end-2) + penalty;
                          cost_p(2:end-1);
                          cost_p(3:end) + penalty;
                          cost_p(4:end) + 4 * penalty diff;
                          cost_p(5:end) + 5 * penalty diff diff], [], 1);
            cost_p = [diff cost(col,2:end-1)+z diff];
            map(col, 2:end-1) = (2:size(cost, 2) - 1) + (idx - 4);
        end

        [~, idx] = min(cost_p);
        disparity_map(row, 1) = idx;

        for col = 1:(cols-1)
            disparity_map(row, col+1) = map(col, max(1, min(size(map,2), round(disparity_map(row,col)))));
        end
    end   
end