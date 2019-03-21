%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outlier_map = calculate_outlier(left_right, right_left, LR_threshold)
    % Initialize outlier map
    [rows, cols] = size(left_right);
    outlier_map = zeros(rows, cols);
    
    for row = 1:rows
        for col = 1:cols
            
            dLR = left_right(row, col);
            
            if col + dLR < 1 || col + dLR > cols                
                outlier_map(row, col) = 1;
            else
                xdLR = round(col + dLR);
                dRL =  right_left(row, xdLR);
                
                if abs(dLR - dRL) <= LR_threshold
                    outlier_map(row,col) = 0;
                else
                    outlier_map(row,col) = 1;
                end
            end
            
        end
    end
end