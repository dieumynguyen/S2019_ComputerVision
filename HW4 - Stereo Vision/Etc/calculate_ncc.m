%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map = calculate_ncc(left_img, right_img, window_size, max_disparity)
    
    % Initialize disparity map & get size of left img
    rows = size(left_img, 1);
    cols = size(left_img, 2);
    
    disparity_map = zeros(rows, cols);

    % Pad right img & im2double
    padding = NaN(rows, max_disparity);
    right_img = horzcat(im2double(right_img), padding);
    left_img = im2double(left_img);

    % Center of window
    center = (window_size - 1) / 2;
    
    % Loop over img rows and cols
    for row = 1 + center:rows - center
        for col = 1 + center:cols - center 
            pix_left = left_img(row-center:row+center, col-center:col+center);
         
            mean_left = mean2(pix_left);
            norm_left = pix_left - mean_left;
            error_left = sqrt(sum(sum((pix_left - mean_left).^2)));
            
            % SD = 0 breaks  NCC
            if error_left == 0
                error_left = 1;
            end
            
            pix_left = repmat(norm_left/error_left, 1, max_disparity);
            
            pix_right = right_img(row-center:row+center, col-center:col+center);
            mean_right = mean2(pix_right);
            norm_right = pix_right - mean_right;
            
            SD_error_right = sqrt(sum(sum((pix_right - mean_left).^2)));
            
            % SD = 0 breaks  NCC
            if SD_error_right == 0
                SD_error_right = 1;
            end
            
            pix_right = norm_right/SD_error_right;

            for d = 1:max_disparity
               c = col + d;
               window = right_img(row-center:row+center, c-center:c+center);
               mean_window = mean2(window);
               norm_mean_window = window-mean_window;
               error_window = sqrt(sum(sum((window-mean_window).^2)));
               
               % SD = 0 breaks  NCC
               if error_window == 0
                   error_window = 1;
               end
               
               window = norm_mean_window/error_window;
               
               pix_right = horzcat(pix_right, window);
            end
            % Make vectors of windows and sum NCC score. 
            NCC = sum(reshape(pix_left.*pix_right, window_size.^2, max_disparity+1));
            
            % Maximum score. 
            disp = find(NCC == max(NCC));
            
            % Current rule: more than one NCC score, taking the first one.
            bestMatch = disp(1);
            
            % Update the disparity map with calculate disparity. 
            disparity_map(row, col) = bestMatch-1;
        end
    end  
end