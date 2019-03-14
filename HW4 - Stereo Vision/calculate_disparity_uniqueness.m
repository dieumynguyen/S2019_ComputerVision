function [disparityMap] = calculate_disparity_uniqueness(leftImage, rightImage, dMin, dMax, window_size, uniqueness_factor)

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
           
           % determine the next lowest SSD score (next best match)
           sorted_SSD = sort(SSD(:));
           disp_2 = sorted_SSD(2);
           
           % Uniqueness constraint: For a match to be unique,
           % the matching cost (SSD) times a uniqueness factor 
           % must be smaller than the cost for the next best match
           % This is like how much the best/min cost value should
           % "win" over the second best to be considered the correct
           % match
           if min(SSD) * uniqueness_factor < disp_2
                bestMatch = disp(1);
                % update disparityMap with calculate disparity. 
                disparityMap(row,col) = bestMatch-1;
           else
                disparityMap = disparityMap;
           end
          
%            % update disparityMap with calculate disparity. 
%            disparityMap(row,col) = bestMatch-1;
       end
   end
end