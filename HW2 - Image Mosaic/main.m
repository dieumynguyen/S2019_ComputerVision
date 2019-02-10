
clear all; close all; clc;

% Random seed
rng(42);

% Read in images 
img_1 = imread('Materials/Square0.jpg');
img_2 = imread('Materials/Square1.jpg');

% Global variables
num_points = 10;

%-------- Task 1. Getting correspondences --------%
% x_1 = [];
% y_1 = [];
% x_2 = [];
% y_2 = [];
% for i = 1:num_points
%    % Select a point from image 1
%    figure(i), imshow(img_1);
%    hold on;
%    [x1, y1] = ginput(1);
%    x_1(i) = x1;
%    y_1(i) = y1;
%    hold off;
%    close(figure(i));
%    
%    % Select the corresponding point in image 2
%    figure(i), imshow(img_2);
%    hold on;
%    [x2, y2] = ginput(1);
%    x_2(i) = x2;
%    y_2(i) = y2;
%    hold off;
%    close(figure(i));
% end

% For debugging
x_1 = [311,262,3.660000000000001e+02,285,3.580000000000001e+02,395,548,531,622,4.610000000000001e+02];
y_1 = [81.999999999999890,1.169999999999999e+02,1.189999999999999e+02,2.949999999999999e+02,2.909999999999999e+02,1.339999999999999e+02,1.789999999999999e+02,1.579999999999999e+02,1.899999999999999e+02,2.839999999999999e+02];
x_2 = [1.060000000000000e+02,55.999999999999960,163,74.999999999999960,148,192,3.390000000000001e+02,325,4.080000000000001e+02,255];
y_2 = [76.999999999999890,1.139999999999999e+02,1.199999999999999e+02,2.969999999999999e+02,2.949999999999999e+02,1.349999999999999e+02,1.829999999999999e+02,1.639999999999999e+02,1.969999999999999e+02,2.849999999999999e+02];

coordinate_points = [x_1', y_1', x_2', y_2'];

%-------- Task 2. Computing the homography parameters --------%
% Ah = b

% Compute A
% Enforcing 8 degree of freedoms
num_H_points = 4;
% A = zeros(2*num_H_points, 8);
% for i = 1:num_H_points
%     A(2*i-1, :) = [x_1(i), y_1(i), 1, 0, 0, 0, -x_1(i)*x_2(i) -y_1(i)*x_2(i)];
%     A(2*i,   :) = [0, 0, 0, x_1(i), y_1(i), 1, x_1(i)*y_2(i) y_1(i)*y_2(i)];
% end
% 
% % Compute b
% b = zeros(2*num_H_points, 1);
% for i = 1:num_H_points
%    b(2*i-1) = x_2(i);
%    b(2*i) = y_2(i);
% end
% 
% % Compute h
% % Source: http://www.cse.psu.edu/~rtc12/CSE486/lecture16.pdf
% h = A\b;
% h(9) = 1;
% H = transpose(reshape(h, [3, 3]));

H = [];
for i = 1:20
    % RANSAC algorithm
    A = zeros(2*num_H_points, 9);

    % Matrix data
    random_index = randsample(1:length(coordinate_points), 4, false);
    random_points = coordinate_points(random_index,:);

    for j = 1 : 4

        x1 = random_points(j, 1);
        y1 = random_points(j, 2);
        x2 = random_points(j, 3);
        y2 = random_points(j, 4);

        A(2*j-1, :) = [x1 y1 1 0 0 0 -x1*x2 -y1*x2 -x2];
        A(2*j,   :) = [0 0 0 x1 y1 1 -x1*y2 -y1*y2 -y2];
    end
    [U,S,V] = svd(A);
    tmp_H = reshape(V(:,end),[3,3]).';
    H = tmp_H / tmp_H(3,3);

    % RANSAC
    % calculate the distance of this H
    min_dis = Inf;
    tmp_dis = 0;
    for j = 1 : 10
        p = [coordinate_points(j,1:2), 1];
        cor_p = [coordinate_points(j,3:4), 1].';
        pro_p = tmp_H * p.';
        pro_p = pro_p / pro_p(3,1);
        tmp_dis = tmp_dis + sum((pro_p-cor_p).^2).^0.5;
    end
    % choose the H with smallest distance
    if tmp_dis < min_dis
        min_dis = tmp_dis;
        H = tmp_H;
    end
end


%-------- Task 3. Warping between image planes --------%

homography = inv(H);

img_1 = im2double(img_1);
img_2 = im2double(img_2);

[X, Y, C] = size(img_2);

% Set the bounding box
bb = [1 Y 1 X];
bb_xmin = bb(1);
bb_xmax = bb(2);
bb_ymin = bb(3);
bb_ymax = bb(4);

[U, V] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
[nrows, ncols] = size(U);

% Compute warped x- and y- coordinates
u = U(:);
v = V(:);
x1 = homography(1, 1)*u + homography(1, 2)*v + homography(1, 3);
y1 = homography(2, 1)*u + homography(2, 2)*v + homography(2, 3);
w1 = 1./(homography(3, 1)*u + homography(3, 2)*v + homography(3, 3));
U(:) = x1 .* w1;
V(:) = y1 .* w1;

% Compute interpolation and set NaN to 0 (black)
warped_img = zeros(nrows, ncols, 3);
warped_img(1:nrows,1:ncols,1) = interp2(img_2(:,:,1), U, V, 'bilinear');
warped_img(1:nrows,1:ncols,2) = interp2(img_2(:,:,2), U, V, 'bilinear');
warped_img(1:nrows,1:ncols,3) = interp2(img_2(:,:,3), U, V, 'bilinear');
warped_img(isnan(warped_img)) = 0;

%-------- Task 4. Create the output mosaic --------%
% Create new image large enough to hold both views

% Overlay one view onto the other
% Leave black where no data is available




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% REDO %%%%%%%%%%%%%%%%%%%%%%%%%
% % Random seed
% rng(42);
% 
% % Read in 2 images
% image1 = imread('Materials/Image1.jpg');
% image2 = imread('Materials/Image2.jpg');
% num_points = 10;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% Task 1: Getting correspondences
% % figure(1);
% % imshow(image1);
% % hAxes1 = gca;
% % title('Image 1');
% % 
% % figure(2);
% % imshow(image2);
% % hAxes2 = gca;
% % title('Image 2');
% % 
% % x1_array = [];
% % y1_array = [];
% % x2_array = [];
% % y2_array = [];
% % count = 1;
% % while (count <= num_points)
% %     axes(hAxes1);
% %     [x1,y1] = ginput(1);
% %     
% %     axes(hAxes2);
% %     [x2,y2] = ginput(1);
% %     count = count + 1;
% %     
% %     x1_array(count) = x1;
% %     y1_array(count) = y1;
% %     x2_array(count) = x2;
% %     y2_array(count) = y2;
% % end
% % 
% 
% % % Get 4 points
% % imshow(image1);
% % [y1, x1] = ginput(4);
% % % % hold on; % Prevent image from being blown away.
% % % % plot(x1, y1, 'r+', 'MarkerSize', 50);
% % % 
% % imshow(image2);
% % [y2, x2] = ginput(4);
% 
% % point_matrix = [x1, y1, x2, y2];
% % point_matrix(1, :) = [];
% % save('point_matrix.mat', 'point_matrix')
% load('point_matrix.mat')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% Task 2: Computing the homography parameters
% 
% % Looping x times to find smallest distance error
% n_trials = 10;
% % Initialize smallest distance to compare iteratively 
% smallest_distance = 1e6;
% % Best H matrix
% best_H = [];
% for i = 1:n_trials
%     % 1. Pick 4 random points from point_matrix
%     random_points = point_matrix(randi([1 4], [1 4]),:,:);
%     num_chosen_points = size(random_points, 1);
% 
%     % 2. Compute 3x3 homographic matrix H
%     A = zeros(2*num_chosen_points, 9);
%     for i = 1:num_chosen_points
%         % Grab a point from matrix
%         point = random_points(i, :);
%         % Image 1 xy coords
%         x1 = point(1);
%         y1 = point(2);
%         % Image 2 xy coords
%         x2 = point(3);
%         y2 = point(4);
% 
%         % Fill in 2 rows in A matrix for this point
%         A(2*i-1,:) = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2, -x2];
%         A(2*i  ,:) = [0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2, -y2];
%     end
% 
%     % Calculate SVD value and obtain 3x3 H matrix
%     [~, ~, V] = svd(A);
%     V_vec = V(:,end);   % Right singular vectors are columns in V
%     H = reshape(V_vec, 3, 3);
% 
%     % 3. Calculate Euclidian distance bw selected points in 1 image,
%     % and the projection of their corresponding points using H
%     image1_coords = random_points(:,[1,2]);
%     image2_coords = random_points(:,[3,4]);
%     third_col = ones(4, 1);
%     image2_coords_homogeneous = [image2_coords third_col];
% 
%     % Find q (corresponding point of selected point)
%     % p is in homogenous coords
%     % q is not. Divide by 3rd value.
%     % q = p * H;
%     p = image2_coords_homogeneous;
%     q = p * H;
%     q_normalization = q(:,3);
%     new_xy = q ./ q_normalization;
%     new_xy(:,3) = [];
% 
%     % Euclidian distance
%     % 4. Save distance value and repeat 20x
%     euclidian_dist = sqrt((image1_coords(:,1) - new_xy(:,1)).^2 + ...
%                           (image1_coords(:,2) - new_xy(:,2)).^2);
%     sum_dist = sum(euclidian_dist);
%     
%     % Comparing smallest distance
%     % Return H matrix with smallest distance error
%     if sum_dist < smallest_distance
%        smallest_distance = sum_dist; 
%        best_H = H;
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% Task 3: Warping between image planes
% % Take the recovered homography matrix & 1 image
% % Return a new image, a warp of the input using H (or its inverse)
% 
% [X, Y, C] = size(image2);
% 
% % Set the bounding box
% bb = [1 Y 1 X];
% bb_xmin = bb(1);
% bb_xmax = bb(2);
% bb_ymin = bb(3);
% bb_ymax = bb(4);
% 
% [U, V] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
% [nrows, ncols] = size(U);
% 
% % Compute warped x- and y- coordinates
% invH = inv(H);
% 
% % u = U(:);
% % v = V(:);
% x1 = invH(1, 1)*U + invH(1, 2)*V + invH(1, 3) ./(invH(3, 1)*U + invH(3, 2)*V + invH(3, 3));
% y1 = invH(2, 1)*U + invH(2, 2)*V + invH(2, 3) ./(invH(3, 1)*U + invH(3, 2)*V + invH(3, 3));
% 
% % Compute interpolation and set NaN to 0 (black)
% image2 = im2double(image2);
% % warped_img(nrows, ncols, 3) = 1;
% % warped_img = zeros(nrows, ncols, 3);
% warped_img(:,:,1) = interp2(image2(:,:,1), x1, y1, 'bilinear');
% warped_img(:,:,2) = interp2(image2(:,:,2), x1, y1, 'bilinear');
% warped_img(:,:,3) = interp2(image2(:,:,3), x1, y1, 'bilinear');
% warped_img(isnan(warped_img)) = 0;
