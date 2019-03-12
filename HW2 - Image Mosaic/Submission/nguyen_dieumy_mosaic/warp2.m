%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [warped_img, up_offset, left_offset] = warp2(img, H)
    % This function was my first attempt to warp
    % Not vectorized or optimal but works ok

    [m, n, ~] = size(img);
    invH = inv(H);

    % Get bounds of the warp
    up_offset = 1e10;
    down_offset = -1e10;
    left_offset = 1e10;
    right_offset = -1e10;
    for i = 1:m
        for j = 1:n
            coords = H * [i,j,1].';
            coords = coords / coords(3,1);
            if coords(1) < up_offset
                up_offset = round(coords(1));
            end
            if coords(1) > down_offset
                down_offset = round(coords(1));
            end
            if coords(2) < left_offset
                left_offset = round(coords(2));
            end
            if coords(2) > right_offset
                right_offset = round(coords(2));
            end
        end
    end

    % Apply invH * (x,y), create warped image
    warped_img = uint8(zeros(down_offset-up_offset, right_offset-left_offset, 3));
    [mm, nn, ~] = size(warped_img);
    for i = 1:mm
        for j = 1:nn
            x = i + up_offset;
            y = j + left_offset;
            coord = invH * [x,y,1].';
            coord = round(coord / coord(3,1));
            if coord(1,1) > 0 && coord(1,1) <= m && coord(2,1) > 0 && coord(2,1) <= n
                warped_img(i,j,:) = img(coord(1,1), coord(2,1),:);
            end
        end
    end
end
