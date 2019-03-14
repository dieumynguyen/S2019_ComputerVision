% Computing the disparity map between 2 images 
% using dynamic programming with backtracking
function I_disp = stereo_dp(left_I, right_I, disp_range, h_block_size, penalty)

    left_I = mean(left_I, 3);
    right_I = mean(right_I, 3);
    
    % Initialize disparity map
    I_disp = zeros(size(left_I), 'single');

    row = size(left_I, 1);
    col = size(left_I, 2);

    % Keep track of cost
    min_diff = Inf;
    disp_cost = min_diff * ones(col, 2*disp_range+1, 'single');

    % Compute disparity using dynamic programming with backtracking
    % Loop over rows
    for m=1:row
        disp_cost(:) = min_diff;
        row_min = max(1, m - h_block_size);
        row_max = min(row, m + h_block_size);

        % Loop over columns
        for n=1:col
            col_min = max(1, n - h_block_size);
            col_max = min(col, n + h_block_size);

            % Set the pixel search limit with disparity range
            pix_min = max(-disp_range, 1 - col_min);
            pix_max = min(disp_range, col - col_max);

            % Search within that range, calc costs
            for i=pix_min:pix_max
                left_pix = left_I(row_min:row_max, (col_min:col_max)+i);
                right_pix = right_I(row_min:row_max, col_min:col_max);
                % Squared error measure bw pixels
                disp_cost(n, i+disp_range+1) = sumsqr(left_pix - right_pix);
            end
        end

        % Compute minimal cost alignment of 2 scanlines 
        % with backtracking 
        Index = zeros(size(disp_cost), 'single');
        c_p = disp_cost(end, :);
        for j=col-1:-1:1
            diff = (col - j + 1) * min_diff;
            [z, inx] = min([diff diff c_p(1:end-4) + 3 * penalty;
                          diff c_p(1:end-3) + 2 * penalty;
                          c_p(1:end-2) + penalty;
                          c_p(2:end-1);
                          c_p(3:end) + penalty;
                          c_p(4:end) + 4 * penalty diff;
                          c_p(5:end) + 5 * penalty diff diff], [], 1);
            c_p = [diff disp_cost(j,2:end-1)+z diff];
            Index(j, 2:end-1) = (2:size(disp_cost, 2) - 1) + (inx - 4);
        end

        [~, inx] = min(c_p);
        I_disp(m, 1) = inx;

        for k = 1:(col-1)
            I_disp(m, k+1) = Index(k, max(1, min(size(Index,2), round(I_disp(m,k)))));
        end

    end
    
end