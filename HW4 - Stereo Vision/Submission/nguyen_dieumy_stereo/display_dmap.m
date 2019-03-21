%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d_color = display_dmap(disparity_map)
    max_d = max(max(disparity_map));
    min_d = min(min(disparity_map));
    
    disparity_map = disparity_map - min_d;
    diff = max_d - min_d;
    
    disparity_map = disparity_map / diff;
    d_color = repmat(disparity_map, 1, 1, 3);
    
    [rows, cols] = size(d_color);
    
    for rpw = 1:rows
        for col = 1:cols
            if isnan(d_color(rpw,col,:))
                d_color(rpw,col,1) = 1;
                d_color(rpw,col,2) = 0; 
                d_color(rpw,col,3) = 0;
            end     
        end
    end

end