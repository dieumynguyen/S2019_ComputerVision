
% Random seed
rng(42);

% Read in 2 images
image1 = imread('Materials/Image1.jpg');
image2 = imread('Materials/Image2.jpg');
num_points = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Task 1: Getting correspondences
% figure(1);
% imshow(image1);
% hAxes1 = gca;
% title('Image 1');
% 
% figure(2);
% imshow(image2);
% hAxes2 = gca;
% title('Image 2');
% 
% x1_array = [];
% y1_array = [];
% x2_array = [];
% y2_array = [];
% count = 1;
% while (count <= num_points)
%     axes(hAxes1);
%     [x1,y1] = ginput(1);
%     
%     axes(hAxes2);
%     [x2,y2] = ginput(1);
%     count = count + 1;
%     
%     x1_array(count) = x1;
%     y1_array(count) = y1;
%     x2_array(count) = x2;
%     y2_array(count) = y2;
% end
% 
% point_matrix = [x1_array' y1_array' x2_array' y2_array'];
% point_matrix(1, :) = [];
% save('point_matrix.mat', 'point_matrix')

load('point_matrix.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Task 2: Computing the homography parameters
% 1. Pick 4 random points from point_matrix
random_points = point_matrix(randi([1 10], [1 4]),:,:);
num_chosen_points = size(random_points, 1);

% 2. Compute 3x3 homographic matrix H
A = zeros(2*num_chosen_points, 9);
for i = 1:num_chosen_points
    % Grab a point from matrix
    point = random_points(i, :);
    % Image 1 xy coords
    x1 = point(1);
    y1 = point(2);
    % Image 2 xy coords
    x2 = point(3);
    y2 = point(4);
    
    % Fill in 2 rows in A matrix for this point
    A(2*i-1,:) = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2, -x2];
    A(2*i  ,:) = [0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2, -y2];
end

% Calculate SVD value and obtain 3x3 H matrix
[~, ~, V] = svd(A);
V_vec = V(:,end);   % Right singular vectors are columns in V
H = reshape(V_vec, 3, 3);

% 3. Calculate Euclidian distance bw selected points in 1 image,
% and the projection of their corresponding points using H
image1_coords = random_points(:,[1,2]);
image2_coords = random_points(:,[3,4]);
third_col = ones(4, 1);
image2_coords_homogeneous = [image2_coords third_col];

% Find q (corresponding point of selected point)
% p is in homogenous coords
% q is not. Divide by 3rd value.
% q = p * H;
p = image2_coords_homogeneous;
q = p * H;
q_normalization = q(:,3);
new_xy = q ./ q_normalization;
new_xy(:,3) = [];

% Euclidian distance
euclidian_dist = sqrt((image1_coords(:,1) - new_xy(:,1)).^2 + ...
                      (image1_coords(:,2) - new_xy(:,2)).^2);

% 4. Save distance value and repeat 20x. 

% Return H matrix with smallest distance error.
