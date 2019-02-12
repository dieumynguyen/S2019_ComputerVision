%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722 - Distance
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = computeH(coordinates_matrix, num_points, num_trials)
    % Initialize some large distance
    min_distance = 1e10;

    % Loop through trials to select best H with lowest distance
    for i = 1:num_trials
        % Select 4 random points
        four_points = datasample(coordinates_matrix, 4, 'Replace', false);

        % Compute A
        A = [];
        for j = 1 : 4
            x1 = four_points(j,1);
            y1 = four_points(j,2);
            x2 = four_points(j,3);
            y2 = four_points(j,4);
            A_j = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2, -x2;
                   0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2, -y2];
            A = cat(1, A, A_j);
        end

        % Compute H
        [~, ~, V] = svd(A);
        H_V = reshape(V(:,end), [3,3]).';
        H_V = H_V / H_V(3,3);

        % Compute Euclidian distance of this H_V
        distance = 0;
        for j = 1:num_points
            % Corresponding points
            p = [coordinates_matrix(j, 1:2), 1];
            corresponding_p = [coordinates_matrix(j, 3:4), 1].';
            % Projection of corresponding points
            projection_p = H_V * p.';
            % Divide by 3rd value to get new xy coords
            projection_p = projection_p / projection_p(3, 1);
            % Compute distance
            distance = distance + sum((projection_p-corresponding_p).^2).^0.5;
        end

        % Select H with smallest distance
        if distance < min_distance
            min_distance = distance;
            H = H_V;
        end
    end
end
