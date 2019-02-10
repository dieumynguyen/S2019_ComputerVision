
% Random seed
rng(42);

% Read in 2 images
image1 = imread('Materials/Image1.jpg');
image2 = imread('Materials/Image2.jpg');
image2_double = im2double(image2);
num_points = 10;

% -------------------------------------
% 1. Select points
num_points = 10;
% for i = 1:num_points
%     figure(i), imshow(image1);
%     hold on;
%     [x, y] = ginput(1);
%     coordinate_points(i, 1) = x;
%     coordinate_points(i, 2) = y;
%     hold off;
%     close(figure(i));
%     
%     figure(i), imshow(image2);
%     hold on;
%     [x, y] = ginput(1);
%     coordinate_points(i, 3) = x;
%     coordinate_points(i, 4) = y;
%     hold off;
%     close(figure(i));
% end

% save('coordinate_points.mat', 'coordinate_points')
load('coordinate_points.mat')

% -------------------------------------
% 2. Compute H & run RANSAC algorithm 
trials = 20;
num_random_points = 4;
minimum_distance = 1e6;
Homography = zeros(3, 3);

for iteration = 1:trials
    random_index = randsample(1:length(coordinate_points), 4, false);
    random_points = coordinate_points(random_index,:);
    
    image1_coords = random_points(:, [1,2]);
    image2_coords = random_points(:, [3,4]);
    
    A = zeros(2*num_random_points, 9);
    
    for i = 1:num_random_points
        x1 = image1_coords(i, 1);
        y1 = image1_coords(i, 2);
        x2 = image2_coords(i, 1);
        y2 = image2_coords(i, 2);

        % Fill in 2 rows in A matrix for this point
        A(2*i-1,:) = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2, -x2];
        A(2*i  ,:) = [0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2, -y2];
    end
    
    % Calculate SVD value and obtain 3x3 H matrix
    [~, ~, V] = svd(A);
    V_vec = V(:,end);   % Right singular vectors are columns in V
    H = reshape(V_vec, 3, 3);

    % Calculate projection points of image 2
    image_coords = coordinate_points(:, [3,4]);
    coords1 = image_coords';
    q = H * [coords1; ones(1, size(coords1,2))];
    p = q(3,:);
    % Normalize back to x,y coords
    projection_point = [q(1,:)./p; q(2,:)./p]; 
    
    % Calculate Euclidian distance
    reference_coords = coordinate_points(:,[1,2]);
    projection_point = projection_point.';
    distance = 0;
    for row = 1 : size(projection_point,1)
        x = reference_coords(row,:);
        c = projection_point(row,:);
        distance = distance + sqrt((x(1) - c(1)).^2 + (x(2) - c(2)).^2);
    end 
    
    % Save H with minimum distance
    if distance < minimum_distance
       minimum_distance = distance;
       Homography = H;
    end
end

% -------------------------------------
% 3. Warp one image
[X, Y, C] = size(image2_double);

% Set the bounding box
bb = [1 Y 1 X];
bb_xmin = bb(1);
bb_xmax = bb(2);
bb_ymin = bb(3);
bb_ymax = bb(4);

[xi, yi] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);

h = inv(Homography);

% simple application of H * (x,y) while changing (x,y) coords to
% homogenous and converting back to (x,y) after translation.
xx = (h(1,1)*xi+h(1,2)*yi+h(1,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
yy = (h(2,1)*xi+h(2,2)*yi+h(2,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));

% solve for warpedImage pixels using interp2. Do each color channel
% separately for RGB images. 
warpedImg(:,:,1) = (interp2(image2_double(:,:,1),xx,yy, 'bilinear',0));
warpedImg(:,:,2) = (interp2(image2_double(:,:,2),xx,yy, 'bilinear',0));
warpedImg(:,:,3) = (interp2(image2_double(:,:,3),xx,yy, 'bilinear',0));

% % Compute warped x- and y- coordinates
% invH = inv(Homography);
% x1 = invH(1, 1)*U + invH(1, 2)*V + invH(1, 3) ./ (invH(3, 1)*U + invH(3, 2)*V + invH(3, 3));
% y1 = invH(2, 1)*U + invH(2, 2)*V + invH(2, 3) ./ (invH(3, 1)*U + invH(3, 2)*V + invH(3, 3));
% 
% % Compute interpolation and set NaN to 0 (black)
% % warped_img = zeros(nrows, ncols, 3);
% warped_img(:,:,1) = interp2(image2_double(:,:,1), U, V, 'bilinear');
% warped_img(:,:,2) = interp2(image2_double(:,:,2), U, V, 'bilinear');
% warped_img(:,:,3) = interp2(image2_double(:,:,3), U, V, 'bilinear');
% % warped_img(isnan(warped_img)) = 0;

