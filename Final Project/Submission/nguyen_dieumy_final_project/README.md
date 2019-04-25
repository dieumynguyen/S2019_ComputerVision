Name: Dieu My Nguyen <br>
Semester: Spring 2019 <br>
Course Number: CSCI 5722 - Distance <br>
Final Project: Detection of honey bees in dense environment <br>
Instructor: Ioana Fleming

The goal of this project is to apply computer vision techniques to detect honey bees from my research data. The data is videos that need preprocessing. For each frame of a video, a image processing will be applied to crop unnecessary parts of the frame, correct uneven illumination, and detect single bees as best as possible. The mesh cage on the top right is where the queen is. Once bees are detected, we will calculate the average distance of all bees to the queen (center of cage). The eventual goal (past the scope of this class project) would be sufficient detection to perform tracking of individual bees over frames to gain some biological insights.

This project is done in Jupyter notebook running Python. The Main.ipynb notebook is the main program that runs the image loading, processing, and writing pipeline, which is encapsulated in these Python files in the modules/ directory:
- config.py: Place to specify data directory and pipeline and visualization parameters (in root directory).
- utils.py: Contains the class VideoHandler to load videos, and the imgs2vid() function to stitch images into a video.
- visualization.py: Contains functions to create visualization of bee detections.
- image_processing.py: Contains functions to process images and perform bee detection.
- pipeline.py: Contains functions to run through a video and apply image processing methods to clean up images and detect bees in them. Also contains code to write a video of processed images.

To run: Jupyter notebook and Python 3.7 are required. Then, simply run the cells in the Main.ipynb notebook. Sample outputs are retained to show examples.

Results/ directory contains results of various free parameters.
