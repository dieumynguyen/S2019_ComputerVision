function [points3D] = reconstruct(disparityMap, stereoParams)
    [rows, cols] = size(disparityMap);
    % create output matrix. 
    points3D = zeros(rows, cols);
    pp1 = stereoParams.CameraParameters1.PrincipalPoint;pp2 = stereoParams.CameraParameters2.PrincipalPoint; 
    % obtain focal length of the left camera 
    F_len = stereoParams.CameraParameters1.FocalLength;
    F_avg = (F_len(1) + F_len(2))/2;
    baseline_distance = sqrt((pp2(1)-pp1(1)).^2 + (pp2(2)-pp1(2)).^2);
    for row = 1 : rows
        for col=1 : cols
            value = ((F_avg * baseline_distance)/disparityMap(row,col))/1000;
            if value >= 1
                value = 1;
            else
                value = 0;
            end
            points3D(row, col) = value;
        end
    end
end