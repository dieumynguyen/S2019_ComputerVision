%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mosaic_img = mosaic(unwarped_img, warped_img, bb_xmin, bb_ymin, f1, f2)

    [m1, n1, ~] = size(unwarped_img);
    [m2, n2, ~] = size(warped_img);

    % Find size of mosaic space using warped and unwarped images
    m = max(m1, m2+bb_xmin);
    if bb_xmin < 0
        m = m - bb_xmin;
    end

    n = max(n1, n2+bb_ymin);
    if bb_ymin < 0
        n = n - bb_ymin;
    end

    % Place warped image onto mosaic
    mosaic_img = double(zeros([m, n, 3]));
    warped = warped_img(1:m2, 1:n2, :);
    mosaic_img(1:m2, 1:n2, :) = warped;

    % Place unwarped img onto mosaic with appropriate offsets
    img_double = im2double(unwarped_img);
    unwarped = img_double(1:m1, 1:n1, :);
    if bb_xmin < 0 && bb_ymin < 0
        mosaic_img(1-bb_xmin+f1:m1-bb_xmin+f1, 1-bb_ymin-f2:n1-bb_ymin-f2, :) = unwarped;
    elseif bb_xmin < 0 && bb_ymin > 0
        mosaic_img(1-bb_xmin+f1:m1-bb_xmin+f1, 1+bb_ymin-f2:n1+bb_ymin-f2, :) = unwarped;
    elseif bb_xmin > 0 && bb_ymin < 0
        mosaic_img(1+bb_xmin+f1:m1+bb_xmin+f1, 1-bb_ymin-f2:n1-bb_ymin-f2, :) = unwarped;
    else
        mosaic_img(1+bb_xmin:m1+bb_xmin, 1+bb_ymin:n1+bb_ymin, :) = unwarped;
    end
end
