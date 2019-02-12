%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mosaic_img = mosaic2(unwarped_image, warped_image, up_offset, left_offset)
    [m1, n1, ~] = size(unwarped_image);
    [m2, n2, ~] = size(warped_image);

    % Determine size of mosaic
    m = max(m1, m2+up_offset);
    if up_offset < 0
        m = m - up_offset;
    end
    n = max(n1,n2+left_offset);
    if left_offset < 0
        n = n - left_offset;
    end

    % Initialize mosaic image
    mosaic_img = uint8(zeros([m,n,3]));

    % Place unwarped image onto mosaic
    for i = 1:m1
        for j = 1:n1
            if up_offset < 0
                x = i - up_offset;
            else
                x = i;
            end
            if left_offset < 0
                y = j - left_offset;
            else
                y = j;
            end
            mosaic_img(x,y,:) = unwarped_image(i,j,:);
        end
    end

    % Place warped image onto mosaic
    for i = 1:m2
        for j = 1:n2
            if up_offset < 0
                x = i;
            else
                x = i + up_offset;
            end
            if left_offset < 0
                y = j;
            else
                y = j + left_offset;
            end
            if sum(warped_image(i,j,:)) ~= 0
                mosaic_img(x,y,:) = warped_image(i,j,:);
            end
        end
    end
end
