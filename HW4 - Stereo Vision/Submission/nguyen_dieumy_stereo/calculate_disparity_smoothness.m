%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map = calculate_disparity_smoothness(left_img, right_img, max_disparity, ... 
                                                        window_size, smoothness_threshold)
                                                    
    % Initialize disparity map & get size of left img as size of map
    disparity_map = zeros(size(left_img), 'single');

    rows = size(left_img, 1);
    cols = size(left_img, 2);

    left_img = mean(left_img, 3);
    right_img = mean(right_img, 3);
    
    d_costs = [];
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
                  disp = sumsqr(template - window);
                  SSDs(idx, 1) = disp; 
              end
              
          d_costs = [d_costs; SSDs];
        end            
    end  
    
    % Loop over each map, compute overall smoothness cost for each
    % Minimize Engergy = SSD + smoothness_cost
    energies = zeros(size(left_img));
    energies_thresh = zeros(size(left_img));
    for d = 1:length(d_costs)
        for row = 1:rows
            for col = 1:cols
                % Compute smoothness term
                % Sum(neighbor pixels p, q) = |d_p - d_q|
                % where p and q are col and col+1 pixels
                smoothness_cost = abs(d(col) - d(col+1));
                energy = d + smoothness_cost;
                energies(d, 1) = energy;
                energies_thresh(d, 1) = d + smoothness_threshold;
            end
        end
    end
    
    % Keep cost + smoothness_cost < cost + smoothness_threshold
    % Lowest = best map
    keep = energies < energies_thresh;
    keep_disps = d_costs(keep);
    [~, costs] = sort(keep_disps);
    match_index = costs(1,1);
    disparity = match_index + pix_min - 1;
    disparity_map(row, col) = disparity;    
end

