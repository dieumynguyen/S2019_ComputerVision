function [disparityMap] = calculate_disparity_ssd(leftImage, rightImage, dMin, dMax, window_size)

    % removing saturation from images. 
    leftImage = im2double(leftImage); rightImage= im2double(rightImage);    

    % Size of the image
    [rows, cols] = size(leftImage);
        
    % zero disparityMap to be the size of leftImage and rightImage
    disparityMap = zeros(rows, cols);
    
    % determine center pixel of window.
    center_pixel = ((window_size-1)/2);
    
    % Image padding to extend the cornors.
    pad = NaN(rows, dMax);
    
    % Right image catenation with Padding
    rightImage = horzcat(rightImage, pad);
    
    % vectorized approach to increase computation speed with window size 1 
    if window_size == 1
         
        for row = 1+center_pixel:rows-center_pixel  
            
            for col = 1+center_pixel:cols-center_pixel
                
                %create disparity vector containing x+disparity 0-64
                right = rightImage(row, col + dMin : col + dMax-1);
                
                left = leftImage(row,col);
                
                % compute squared differences.
                SSD = (left - right).^2;
                
                % locate the lowest SSD score. 
                disp = find(SSD == min(SSD));
                
                % Best match disparty at this location.
                % Rule: should there be more than one index with min(vec2)
                % use first indexed location
                bestMatch = disp(1);
               
                %place in disparity map.
                disparityMap(row,col) = bestMatch-1;
            end
        end
    end
    if window_size > 1
        % loop through image.
       for row = 1 + center_pixel : rows - center_pixel
           for col = 1 + center_pixel : cols - center_pixel
               left = leftImage(row - center_pixel : row + center_pixel, col - center_pixel : col + center_pixel);
               left = repmat(left,1,65);
               right = rightImage(row - center_pixel : row + center_pixel, col - center_pixel : col + center_pixel);
               for d = 1:64
                  c = col+d;
                  window = rightImage(row - center_pixel : row + center_pixel, c - center_pixel : c + center_pixel);
                  right = horzcat(right, window);
               end
               SSD = (left - right).^2;
               SSD = imgaussfilt(SSD, 'FilterSize', window_size);
               SSD = reshape(SSD, window_size.^2, 65);

               % sum disparities by window. 
               SSD = sum(SSD);
               
               % determine lowest SSD score. 
               disp = find(SSD == min(SSD));
               
               % Rule: if more than 1 best matching disparity then choose
               % first indexed disparity. 
               bestMatch = disp(1);
               
               % update disparityMap with calculate disparity. 
               disparityMap(row,col) = bestMatch-1;
           end
       end
    end
end