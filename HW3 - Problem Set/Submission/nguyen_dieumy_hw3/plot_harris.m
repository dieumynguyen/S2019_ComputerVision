%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 3 Question 9
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For opening image, calling harris(), and plotting corners over img
function plot_harris(im, w, threshold, suppression)
    % Convert image to black and white
    bw = double(im(:,:,1)) ./ 256;

    % Call harris corner detector function
    [corner_coords, descriptors] = harris(bw, w, threshold, suppression);

    % Plot result with marked corners
    imshow(bw);
       
    h = gca;
    h.Visible = 'On';
    if suppression == 1
       title('With non-max suppression');
    else
       title('Without non-max suppression');
    end
    hold on;
    plot(corner_coords(:,1), corner_coords(:,2), 'g+', ...
        'MarkerSize', 12,...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r');   

end