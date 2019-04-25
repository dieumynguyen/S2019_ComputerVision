import matplotlib.pyplot as plt
import cv2

DATASET_DIR = "/video_data/"
VID_OUTPUT_DIR = "processed_vids"

VID_START_FRAME_I = 0
VID_STOP_FRAME_I = 10

REMOVE_BACKGROUND_LINES = True

GRAY_CMAP = plt.cm.gray

DEFAULT_PROCESSING_PARAMS = {
    "SINGLE_BEE_AREA_THRESHOLD" : 1200,
    "BEE_AREA_PERCENTILE"       : 83,        # 85 good
    "LOCAL_THRESHOLDING"        : False,
    "GLOBAL_CROP"               : True,
    "RECURSE_THRESHOLDING"      : 'global',  # global, local
    "NUM_RECURSIONS"            : 10,
    "MAX_FACTOR"                : 0.15,      # 0.15
    "GLOBAL_TRHESHOLD"          : 70,        # 70
    "RECURSE_GLOBAL_THRESHOLD"  : 40,        # 40
    "NON_MAX_SUPRESSION"        : False,
    "NMS_THRESHOLD"             : 0.99,
    "MOROPHOLOGY_TRANSFORM"     : "erode"    # cv2 functions: erode, dilate, MORPH_OPEN, MORPH_CLOSE
    # None, opening, binary_erosion, binary_opening, binary_closing
}

# Visualization Params
OVERLAY_PARAMS = {
    "box_size"         : 2,
    "worker_dot_size"  : 5,
    "queen_dot_size"   : 20,
    "line_thickness"   : 1,
    "font"             : cv2.FONT_HERSHEY_SIMPLEX,
    "font_size"        : 1,
    "cmap"             : plt.cm.cool
}

OVERLAY_FLAGS = {
    "show_worker_dots" : True,
    "show_queen_dot"   : True,
    "show_boxes"       : False,
    "show_lines"       : False,
    "show_text"        : False,
    "show_main_text"   : True
}
