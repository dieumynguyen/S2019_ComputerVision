import cv2
import numpy as np
import skimage.morphology as morphology # opening, binary_erosion, binary_opening, binary_closing

def source_img_preprocessing(img_src_BGR, remove_background_lines):
    # Convert to gray-scale
    img_GRAY = cv2.cvtColor(img_src_BGR, cv2.COLOR_BGR2GRAY)

    # Make copy of source image
    if remove_background_lines:
        img_GRAY = remove_lines(img_GRAY)
        BGR_img = remove_lines_BGR(img_src_BGR)

        # Initialize the overlay image
        overlay_img_src = np.copy(BGR_img)
    else:
        overlay_img_src = np.copy(img_src_BGR)

    return img_GRAY, overlay_img_src

def otsu_thresholding(img):
    _, img_otsu = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    # Invert the black and white?
    img_otsu = cv2.bitwise_not(img_otsu)
    return img_otsu

def global_thresholding(img, global_threshold):
    thresh_img = img < global_threshold
    return thresh_img

def remove_lines(gray):
    rowvals = np.mean(np.sort(gray, axis=1)[:,-100:], axis=1)
    graymod = np.copy(gray).astype(np.float)
    graymod *= np.expand_dims(np.max(rowvals) / np.array(rowvals), axis=1)

    return graymod.astype(np.uint8)

def remove_lines_BGR(img_BGR):
    img_BGR_no_lines = np.zeros_like(img_BGR)
    for c_i in range(3):
        img_BGR_no_lines[:,:,c_i] = remove_lines(img_BGR[:,:,c_i])

    return img_BGR_no_lines

def morphology_transform(img, transform_key):

    kernel = np.ones((3, 3), np.uint8)

    img = (img * 1).astype('uint8')

    # For opening and closing
    if transform_key == "MORPH_OPEN" or transform_key == "MORPH_CLOSE":
        transform = cv2.__dict__[transform_key]
        img_transformed = cv2.morphologyEx(img, transform, kernel)

    else:
        # For only erosion and dilation
        transform = cv2.__dict__[transform_key]
        img_transformed = transform(img, kernel, iterations=1)

    return img_transformed

# Malisiewicz et al.
# https://www.pyimagesearch.com/2015/02/16/faster-non-maximum-suppression-python/
def non_max_suppression_fast(boxes, overlapThresh):
    # if there are no boxes, return an empty list
    if len(boxes) == 0:
        return []

    # if the bounding boxes integers, convert them to floats --
    # this is important since we'll be doing a bunch of divisions
    if boxes.dtype.kind == "i":
        boxes = boxes.astype("float")

    # initialize the list of picked indexes
    pick = []

    # grab the coordinates of the bounding boxes
    x1 = boxes[:,0]
    y1 = boxes[:,1]
    x2 = boxes[:,2]
    y2 = boxes[:,3]

    # compute the area of the bounding boxes and sort the bounding
    # boxes by the bottom-right y-coordinate of the bounding box
    area = (x2 - x1 + 1) * (y2 - y1 + 1)
    idxs = np.argsort(y2)

    # keep looping while some indexes still remain in the indexes
    # list
    while len(idxs) > 0:
        # grab the last index in the indexes list and add the
        # index value to the list of picked indexes
        last = len(idxs) - 1
        i = idxs[last]
        pick.append(i)

        # find the largest (x, y) coordinates for the start of
        # the bounding box and the smallest (x, y) coordinates
        # for the end of the bounding box
        xx1 = np.maximum(x1[i], x1[idxs[:last]])
        yy1 = np.maximum(y1[i], y1[idxs[:last]])
        xx2 = np.minimum(x2[i], x2[idxs[:last]])
        yy2 = np.minimum(y2[i], y2[idxs[:last]])

        # compute the width and height of the bounding box
        w = np.maximum(0, xx2 - xx1 + 1)
        h = np.maximum(0, yy2 - yy1 + 1)

        # compute the ratio of overlap
        overlap = (w * h) / area[idxs[:last]]

        # delete all indexes from the index list that have
        idxs = np.delete(idxs, np.concatenate(([last],
            np.where(overlap > overlapThresh)[0])))

    # return only the bounding boxes that were picked using the
    # integer data type
    return boxes[pick].astype("int"), pick
