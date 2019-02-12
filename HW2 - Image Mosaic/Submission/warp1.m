%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Dieu My Nguyen
% Semester: Spring 2019
% Course Number: CSCI 5722B
% Assignment: 2
% Instructor: Ioana Fleming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [warped_img, bb_xmin, bb_ymin] = warp1(img, H)
    % Warp an image into the plane of another image
    % with a recovered homography matrix, H.

    [m, n, ~] = size(img);

    % Make bounding box
    y = H * [[1;1;1], [1;m;1], [n;m;1], [n;1;1]];
    y(1,:) = y(1,:) ./ y(3,:);
    y(2,:) = y(2,:) ./ y(3,:);

    bb = [ceil(min(y(1,:)));
          ceil(max(y(1,:)));
          ceil(min(y(2,:)));
          ceil(max(y(2,:)));
          ];

    bb_xmin = bb(1);
    bb_xmax = bb(2);
    bb_ymin = bb(3);
    bb_ymax = bb(4);

    [U, V] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
    [nrows, ncols] = size(U);

    % Inverse of H for inverse mapping
    invH = inv(H);

    % Change xy coords to homogeneous and convert back
    u = U(:);
    v = V(:);
    x1 = invH(1,1)*u + invH(1,2)*v + invH(1,3);
    y1 = invH(2,1)*u + invH(2,2)*v + invH(2,3);
    w1 = 1 ./ (invH(3,1)*u + invH(3,2)*v + invH(3,3));
    U(:) = x1 .* w1;
    V(:) = y1 .* w1;

    warped_img(nrows, ncols, 3) = 1;
    img_double = im2double(img);

    % Create warped image, each channel separately for RGB
    warped_img(:,:,1) = interp2(img_double(:,:,1),U,V,'bilinear');
    warped_img(:,:,2) = interp2(img_double(:,:,2),U,V,'bilinear');
    warped_img(:,:,3) = interp2(img_double(:,:,3),U,V,'bilinear');
end
