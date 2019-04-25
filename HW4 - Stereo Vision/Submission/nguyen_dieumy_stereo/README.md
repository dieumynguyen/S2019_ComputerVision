Name: Dieu My Nguyen <br>
Semester: Spring 2019 <br>
Course Number: CSCI 5722 - Distance <br>
Assignment 4: Stereo and Segmentation <br>
Instructor: Ioana Fleming

Task 1 to Task 6 are implemented inside the provided file DepthEstimationFromStereoVideoExample.m using the same images provided it in. Other parts of this file are commented out to only run the needed parts. Task 7 tests the code on the provided on the frame_1L.png and frame_1R.png images. The dynamic programming code is also tested on the frame_1L.png and frame_1R.png images.

<b> Task 1: </b> Calculate disparity using the SSD algorithm. <br>
<b> Task 2: </b> Calculate disparity using the NCC algorithm. <br>
<b> Task 3: </b> Uniqueness constraint. For a stereo match to be considered unique, the minimum matching cost times a uniqueness factor must be smaller than the cost for the next best match. <br>
<b> Task 4: </b> Smoothness constraint. Treat the problem of finding the disparity map as minimizing energy, which is the sum of the data cost and the smoothness cost, defined by the sum of the absolute difference between the disparity of pixel p and of pixel q. <br>
<b> Task 5: </b> Generate outliers map. <br>
<b> Task 6: </b> Compute depth from disparity. <br>
<b> Task 7: </b> Synthetic stereo sequences. <br>
<b> Task 8-10: </b> Dynamic programming to calculate disparity map. <br>


