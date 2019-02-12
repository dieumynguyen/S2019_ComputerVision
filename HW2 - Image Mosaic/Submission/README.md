Name: Dieu My Nguyen
Semester: Spring 2019
Course Number: CSCI 5722 - Distance
Assignment 2: Image Mosaics
Instructor: Ioana Fleming

This is an implementation of an image stitcher that uses image warping and homographies to create an image mosaic from two input images. This is composed of the following functions in their respective scripts:

- getPoint(): Takes 2 input images. Uses the ginput function to manually collect mouse click points that appear in both input images. Outputs a 10x4 matrix, where each row is a pair of corresponding
points, and the 4 columns represent the (row, column) position from the 1st image, followed by the (row, column) position from the 2nd image.

- computeH(): takes the set of corresponding image points
returned by the function getPoints() from Task 1, and computes the associated 3x3 homography matrix H using the RANSAC algorithm. The homography matrix with the smallest distance error is returned.

- warp1(): Takes the recovered homography matrix and one
image, and return a new image that is the (inverse) warp of the input image using H and bilinear interpolation.

- mosaic(): Stitch together the unwarped and warped images by simple overlaying based on matrix indices.

Two preliminary functions, warp2() and mosaic2(), were my first attempt to solve this problem. They're not optimal due to nested for loops, but they also work ok.

Run the complete code in main.m file to obtain 4 figures for the 4 sets of images and their resulting mosaics: the provided Squared images, my own two sets of images (Bee and Hike), and the Frame images. Running main.m will show the 4 figures, each with the 2 original images and their mosaic.

Input images are from Materials/.
The results (4 figures and 4 mosaics) are saved in Results/.  
