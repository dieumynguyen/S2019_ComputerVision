import numpy as np
import cv2
import matplotlib.pyplot as plt

def get_bboxes(frame):
    bboxes = frame['stats']['global']
    new_bboxes = []
    for bbox in bboxes:
        x1 = bbox[0]
        x2 = bbox[0] + bbox[2]
        y1 = bbox[1]
        y2 = bbox[1] + bbox[3]
        new_bboxes.append(np.array([x1, x2, y1, y2]))

    return np.array(new_bboxes)

def create_overlay_img(img_original, bees, config, processing_params, queen_x=0, queen_y=0, show=False, nms=True, nms_threshold=0.95):

    stats = bees['stats']['global']
    centroids = bees['centroids']
    label_color = bees['label_color']

    # Image checking
    # ============================================================
    # Copy image and set to byte
    overlay_img = np.copy(img_original).astype(np.uint8)

    # Generate overlay image
    if len(img_original.shape) == 2:
        overlay_img = cv2.cvtColor(overlay_img, cv2.COLOR_GRAY2RGB)

    # Ensure image in correct range
    if overlay_img.max() == 1:
        overlay_img *= 255
    # ============================================================

    # Queen
    if config.OVERLAY_FLAGS['show_queen_dot']:
        cv2.circle(overlay_img, (queen_x, queen_y), config.OVERLAY_PARAMS['queen_dot_size'], (0, 0, 255), -1)

    # Workers
    distances = []
    for stat, centroid in zip(stats, centroids):
        area = stat[-1]
        top_left = tuple(stat[:2])
        h,w = stat[2:4]
        bottom_right = (top_left[0]+h, top_left[1]+w)
        centroid_x, centroid_y = centroid.astype(np.int)

        if config.OVERLAY_FLAGS['show_lines']:
            cv2.line(overlay_img, (centroid_x, centroid_y), (queen_x, queen_y), (102, 0, 255), config.OVERLAY_PARAMS['line_thickness'])

        if config.OVERLAY_FLAGS['show_boxes']:
#             if area > processing_params['SINGLE_BEE_AREA_THRESHOLD']:
#                 label_color = (0, 0, 255)
#             elif area > processing_params['SINGLE_BEE_AREA_THRESHOLD'] * 1.2:
#                 label_color = (255, 0, 0)
#             elif area > processing_params['SINGLE_BEE_AREA_THRESHOLD'] * 1.5:
#                 label_color = (255, 255, 0)

            label_color = (255, 0, 0)

            cv2.rectangle(overlay_img, top_left, bottom_right, label_color, config.OVERLAY_PARAMS['box_size'])

        if config.OVERLAY_FLAGS['show_worker_dots']:
            label_color = (0, 255, 0)
            cv2.circle(overlay_img, (centroid_x, centroid_y),
                       config.OVERLAY_PARAMS['worker_dot_size'], label_color, -1)

        dx = queen_x-centroid_x
        dy = queen_y-centroid_y
        distance_i = int(np.sqrt((dx)**2 + (dy)**2))
        distances.append(distance_i)

        if config.OVERLAY_FLAGS['show_text']:
            distance_str = f"D2Queen: {distance_i} px"
            cv2.putText(overlay_img, distance_str, (centroid_x, centroid_y),
                        config.OVERLAY_PARAMS['font'], config.OVERLAY_PARAMS['font_size'], (255, 255, 255), 2, cv2.LINE_AA)


            area_str = f"Area: {area} px"
            cv2.putText(overlay_img, area_str, (centroid_x, centroid_y-30),
                        config.OVERLAY_PARAMS['font'], config.OVERLAY_PARAMS['font_size'], (255, 255, 255), 2, cv2.LINE_AA)

    if config.OVERLAY_FLAGS['show_main_text']:
        # Calculate average distance
        ave_distance = int(np.mean(distances))
        ave_distance_str =   f"Avg Distance:   {ave_distance} px"

        # Num detections
        num_detections_str = f"Num Detections: {len(centroids)}"

        cv2.putText(overlay_img, num_detections_str, (150, 150), config.OVERLAY_PARAMS['font'],
                    2, (255, 255, 255), 5, cv2.LINE_AA)
        cv2.putText(overlay_img, ave_distance_str, (150, 210), config.OVERLAY_PARAMS['font'],
                    2, (255, 255, 255), 5, cv2.LINE_AA)

    if show:
        fig, ax = plt.subplots(figsize=(12,12))
        ax.imshow(overlay_img)

    return overlay_img
