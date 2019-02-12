%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outimg] = mosaic2(img1,img2,up,left)
    [m1,n1,k] = size(img1);
    [m2,n2,k] = size(img2);
    m = max(m1,m2+up);
    if up < 0
        m = m - up;
    end
    n = max(n1,n2+left);
    if left < 0
        n = n - left;
    end
    outimg = uint8(zeros([m,n,3]));
    % first, draw the unwarped image
    for i = 1 : m1
        for j = 1 : n1
            if up < 0
                x = i - up;
            else
                x = i;
            end
            if left < 0
                y = j - left;
            else
                y = j;
            end
            outimg(x,y,:) = img1(i,j,:);
        end
    end
    % then, draw the warped image including the overlapped part
    for i = 1 : m2
        for j = 1 : n2
            if up < 0
                x = i;
            else
                x = i + up;
            end
            if left < 0
                y = j;
            else
                y = j + left;
            end
            if sum(img2(i,j,:)) ~= 0
                outimg(x,y,:) = img2(i,j,:);
            end
        end
    end
end

