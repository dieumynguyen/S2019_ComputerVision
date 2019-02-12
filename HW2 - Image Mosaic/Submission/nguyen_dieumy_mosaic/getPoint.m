%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coordinates_matrix = getPoint(img1, img2, num_points)
    figure
    imagesc(img1);
    [y1, x1] = ginput(num_points);
    imagesc(img2);
    [y2, x2] = ginput(num_points);
    coordinates_matrix = [x1, y1, x2, y2];
end

