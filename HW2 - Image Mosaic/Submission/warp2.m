%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outimg,up,left] = warp2(inimg, H)
    [m,n,k] = size(inimg);
    inv_H = inv(H);
    % find the boundary of the warp of the image
    up = Inf;
    down = -Inf;
    left = Inf;
    right = -Inf;
    for i = 1 : m
        for j = 1 : n
            tmp = H*[i,j,1].';
            tmp = tmp/tmp(3,1);
            if tmp(1) < up
                up = round(tmp(1));
            end
            if tmp(1) > down
                down = round(tmp(1));
            end
            if tmp(2) < left
                left = round(tmp(2));
            end
            if tmp(2) > right
                right = round(tmp(2));
            end
        end
    end

    outimg = uint8(zeros(down-up,right-left,3));
    [mm nn k] = size(outimg);
    for i = 1 : mm
        for j = 1 : nn
            x = i + up;
            y = j + left;
            sample_cor = inv_H * [x,y,1].';
            sample_cor = round(sample_cor / sample_cor(3,1));
            if sample_cor(1,1) > 0 && sample_cor(1,1) <= m && sample_cor(2,1) > 0 && sample_cor(2,1) <= n
                outimg(i,j,:) = inimg(sample_cor(1,1),sample_cor(2,1),:);
            end
        end
    end
end

