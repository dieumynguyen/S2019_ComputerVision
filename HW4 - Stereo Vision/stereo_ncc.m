function disparityMap = stereo_ssd(leftImage, rightImage, dMin, dMax, window_size)
    
    % remove color saturation.
    leftImage = im2double(leftImage);
        
    % determine size for disparity map.
    [rows, cols] = size(leftImage);
    
    % allocate disparityMap
    disparityMap = zeros(rows, cols);
        
    % image padding.
    pad = NaN(rows, dMax);

    rightImage = horzcat(im2double(rightImage), pad);
    
    % determine distance to center of window
    center_pixel = (window_size-1)/2;
    
    for row = 1 + center_pixel : rows - center_pixel
        
        for col = 1 + center_pixel : cols - center_pixel
             
            left = leftImage(row - center_pixel : row + center_pixel, col - center_pixel : col + center_pixel);
         
            % mean Calculation.
            mean_left = mean2(left);
            
            normalization_left = left - mean_left;
            
            SD_error_left = sqrt(sum(sum((left - mean_left).^2)));
            
            % SD = 0 breaks  NCC
            if SD_error_left == 0
                SD_error_left = 1;
            end
            
            left = repmat(normalization_left/SD_error_left, 1, 65);
            
            right = rightImage(row - center_pixel : row + center_pixel, col - center_pixel : col + center_pixel);
            
            mean_right = mean2(right);
            
            Normalization_right = right - mean_right;
            
            SD_error_right = sqrt(sum(sum((right - mean_left).^2)));
            
            % SD = 0 breaks  NCC
            if SD_error_right == 0
                SD_error_right = 1;
            end
            
            right = Normalization_right/SD_error_right;

            for d = 1:64
                
               c = col + d;
               
               window = rightImage(row - center_pixel:row + center_pixel, c - center_pixel : c + center_pixel);
               
               mean_window = mean2(window);
               
               normalization_mean_window = window-mean_window;
               
               SD_error_window = sqrt(sum(sum((window-mean_window).^2)));
               
               % SD = 0 breaks  NCC
               if SD_error_window == 0
                   SD_error_window = 1;
               end
               
               window = normalization_mean_window/SD_error_window;
               
               right = horzcat(right, window);
            end
            
            % reshape to create vectors of windows and sum NCC score. 
            NCC = sum(reshape(left.*right, window_size.^2, 65));
            
            % determine the maximum score. 
            disp = find(NCC == max(NCC));
            
            % Current rule: more than one NCC score, taking the first one.
            bestMatch = disp(1);
            
            % update the disparity map with calculate disparity. 
            disparityMap(row, col) = bestMatch-1;
        end
    end  
end