%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen 
% Semester: Spring 2019 
% Course Number: CSCI 5722 - Distance 
% Assignment 4: Stereo and Segmentation 
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function points3D = compute_depth(disparity_map, stereo_params)
    [rows, cols] = size(disparity_map);
    
    % Initialize & get cam params
    points3D = zeros(rows, cols);
    pp1 = stereo_params.CameraParameters1.PrincipalPoint;pp2 = stereo_params.CameraParameters2.PrincipalPoint; 
    focal_len = stereo_params.CameraParameters1.FocalLength;
    focal_avg = (focal_len(1) + focal_len(2)) / 2;
    baseline = sqrt((pp2(1)-pp1(1)).^2 + (pp2(2)-pp1(2)).^2);
    
    for row = 1:rows
        for col = 1:cols
            value = ((focal_avg * baseline) / disparity_map(row,col)) / 1000;
            if value >= 1, value = 1;else, value = 0; end
            points3D(row, col) = value;
        end
    end
end