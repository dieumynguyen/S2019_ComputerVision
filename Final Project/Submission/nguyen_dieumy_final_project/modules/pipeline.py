import os
import sys
import numpy as np
import cv2
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

import config

# Import Python files of low-level implementations
import modules.image_processing as image_processing
import modules.visualization as visualization
import modules.utils as utils

######## RUN CONNECTED COMPONENTS

def run_connected_components(img, processing_params, prev_stats=None, preview=False):
    img = np.copy(img).astype(np.uint8)

    # Run connected components
    num_regions, regions, stats, centroids = cv2.connectedComponentsWithStats(img)

    areas = stats[:,-1]

    if preview:
        mean = np.mean(areas)
        median = np.median(areas)
        percentile = np.percentile(areas, processing_params['BEE_AREA_PERCENTILE'])
        plt.figure()
        plt.hist(areas, density=True)
        plt.title(f"Mean: {mean} -- Median: {median} -- 50 Percentile: {percentile}")
        plt.figure()

    # Filter out by area
    MIN_AREA = np.percentile(areas, processing_params['BEE_AREA_PERCENTILE'])
    MAX_AREA = np.product(img.shape[:2])*processing_params['MAX_FACTOR']

    # Select elements where: MIN_AREA < element area < MAX_AREA
    condition = (stats[:,-1] < MAX_AREA) & (stats[:,-1] > MIN_AREA)
    new_stats = stats[condition]
    new_centroids = centroids[condition]

    # Catch case where no elements satisfy constraints
    if len(new_centroids) == 0:
        return None, None

    # Check for previous stats; update global positions if exist
    if prev_stats is None:
        stats = {
            "local"  : new_stats,
            "global" : new_stats
        }
        centroids = new_centroids
    else:
        global_stats = np.copy(new_stats)
        for crop_i in range(len(global_stats)):
            global_stats[crop_i,:2] += prev_stats[:2]
        stats = {
            "local"  : new_stats,
            "global" : global_stats
        }
        centroids = new_centroids
        for centroid_i in range(len(centroids)):
            centroids[centroid_i][0] += prev_stats[0]
            centroids[centroid_i][1] += prev_stats[1]

    return stats, centroids

##########################################################################################

######## CROP BEES

def crop_bees(img, stats):
    src_img = np.copy(img)
    cropped_imgs = []
    for stat in stats:
        top_left_x, top_left_y, width, height = stat[:4]
        cropped_img = src_img[top_left_y:top_left_y+height, top_left_x:top_left_x+width]
        cropped_imgs.append(cropped_img)
    try:
        cropped_imgs = np.array(cropped_imgs)
    except:
        cropped_imgs = None

    return cropped_imgs

##########################################################################################

######## SORT AND UPDATE BEES

def sort_bees(stats, centroids, processing_params):
    good_bees = {"stats" : {"global" : [], "local" : []}, "centroids" : []}
    bad_bees = {"stats" : {"global" : [], "local" : []}, "centroids" : []}

    for i in range(len(centroids)):
        centroid = centroids[i]
        local_stat = stats['local'][i]
        global_stat = stats['global'][i]
        area = global_stat[-1]

        if area <= processing_params['SINGLE_BEE_AREA_THRESHOLD']:
            container = good_bees
        else:
            container = bad_bees

        container['stats']['global'].append(global_stat)
        container['stats']['local'].append(local_stat)
        container['centroids'].append(centroid)

    if good_bees['stats']['global'] == []:
        good_bees = None

    if bad_bees['stats']['global'] == []:
        bad_bees = None

    return good_bees, bad_bees

def update_bees(container, bees, processing_params, level_i=0):
    local_bee_stats = np.array(bees['stats']['local'])
    global_bee_stats = np.array(bees['stats']['global'])
    bee_centroids = np.array(bees['centroids'])

    if container['stats']['local'] is not None:
        local_bee_stats = np.concatenate([container['stats']['local'], local_bee_stats], axis=0)
        global_bee_stats = np.concatenate([container['stats']['global'], global_bee_stats], axis=0)
        bee_centroids = np.concatenate([container['centroids'], bee_centroids], axis=0)

    if level_i == 0:
        label_color = (0, 255, 0)
    else:
        cmap = config.OVERLAY_PARAMS['cmap']
        label_color = np.array(cmap(level_i/processing_params['NUM_RECURSIONS'])[:3])*255

    # Update container
    container['stats']['local'] = local_bee_stats
    container['stats']['global'] = global_bee_stats
    container['centroids'] = bee_centroids
    container['label_color'] = label_color

##########################################################################################

######## BOUNDING BOXES FOR DETECTED BEES

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

##########################################################################################

######## FILTER BEES WITH NON-MAX SUPPRESSION

def filter_bees_by_NMS(bees_i, idxs):
    new_bees_i = {}
    for key, val in bees_i.items():
        if key == 'stats':
            new_bees_i[key] = {}
            new_bees_i[key]['global'] = val['global'][idxs]
        elif key == 'centroids':
            new_bees_i[key] = val[idxs]
        else:
            new_bees_i[key] = val
    return new_bees_i

##########################################################################################

######## RUN MAIN PROCESSING FUNCTION

def process_image(img_GRAY, processing_params):
    # All good bees container
    all_good_bees = {"stats" : {"global" : None, "local" : None}, "centroids" : None}

    src_img = np.copy(img_GRAY)
    original_img = np.copy(img_GRAY)

    # Threshold Image
    # -------------------------------------------------
    if processing_params['LOCAL_THRESHOLDING']:
        src_img = image_processing.otsu_thresholding(src_img)
    else:
        init_global_threshold = processing_params['GLOBAL_TRHESHOLD']
        src_img = image_processing.global_thresholding(src_img, global_threshold=init_global_threshold)
    # -------------------------------------------------

    # Run initial connected components
    new_stats, new_centroids = run_connected_components(src_img, processing_params)

    if new_centroids is None:
        return all_good_bees

    # Sort results into good and bad results - bad results are recursed on
    good_bees, bad_bees = sort_bees(new_stats, new_centroids, processing_params)

    # Update good bees container
    # If new good bees, update all good bees container
    if good_bees is not None:
        update_bees(all_good_bees, good_bees, processing_params, level_i=0)

    if bad_bees is None:
        return all_good_bees
    else:
        update_bees(all_good_bees, bad_bees, processing_params, level_i=255)

    # Begin recursion on bad bees
    for i in range(processing_params['NUM_RECURSIONS']):
        # Extract stats and centroids from all bad bees
        new_stats = bad_bees['stats']
        new_centroids = bad_bees['centroids']

        # Crop bad bee regions from src image
        if processing_params['GLOBAL_CROP']:
            cropped_imgs = crop_bees(original_img, new_stats['global'])
        else:
            cropped_imgs = crop_bees(src_img, new_stats['local'])

        # Check for case where no detections pass thresholding reqs in crop
        if cropped_imgs is None:
            continue

        # Init container for bad bees
        new_bad_bee_batch = {"stats" : {"global" : None, "local" : None}, "centroids" : None}
        for crop_i, cropped_img in enumerate(cropped_imgs):

            # Get centroid and stats for current bad bee crop
            prev_centroids = new_centroids[crop_i]
            prev_stats = {
                "local"  : new_stats['local'][crop_i],
                "global" : new_stats['global'][crop_i],
            }

            # Set source image to crop img
            src_img = np.copy(cropped_img)

            # Threshold image
            if processing_params['RECURSE_THRESHOLDING'] == 'local':
                src_img = image_processing.otsu_thresholding(src_img.astype(np.uint8))
            else:
                area_i = prev_stats['local'][-1] # maybe use area to create new global threshold
                new_global_thresh = processing_params['RECURSE_GLOBAL_THRESHOLD'] / (i+1)
                src_img = image_processing.global_thresholding(src_img, global_threshold=new_global_thresh)

            if processing_params['MOROPHOLOGY_TRANSFORM'] is not None:
                src_img = image_processing.morphology_transform(src_img, processing_params['MOROPHOLOGY_TRANSFORM'])

            # Connected components
            new_stats_i, new_centroids_i = run_connected_components(src_img, processing_params, prev_stats=prev_stats['global'])

            # If no centroids are found, append bad bee final result to all good bees
            if new_centroids_i is None:
                prev_stats['local'] = np.expand_dims(prev_stats['local'], axis=0)
                prev_stats['global'] = np.expand_dims(prev_stats['global'], axis=0)
                prev_centroids = np.expand_dims(prev_centroids, axis=0)
                bad_bees_end = {
                    "stats"     : prev_stats,
                    "centroids" : prev_centroids
                }
                update_bees(all_good_bees, bad_bees_end, processing_params, level_i=i+1)
                continue

            # Else, new connected components found
            # Sort bees into new good and new bad
            new_good_bees, new_bad_bees = sort_bees(new_stats_i, new_centroids_i, processing_params)

            # If new good bees, update all good bees container
            if new_good_bees is not None:
                update_bees(all_good_bees, new_good_bees, processing_params, level_i=i+1)

            # If new bad bees, update the batch_i bad bee holder
            if new_bad_bees is not None:
                update_bees(new_bad_bee_batch, new_bad_bees, processing_params, level_i=i+1)

        # Check if any new updates in bad bees - update bad bees list; otherwise, break
        if new_bad_bee_batch['centroids'] is not None:
            bad_bees = new_bad_bee_batch
        else:
            bad_bees = None
            break

    if bad_bees is not None:
        update_bees(all_good_bees, bad_bees, processing_params, level_i=i+1)

    if processing_params["NON_MAX_SUPRESSION"]:
        nms_threshold = processing_params["NMS_THRESHOLD"]
        bboxes = get_bboxes(all_good_bees)
        new_bboxes, idxs = image_processing.non_max_suppression_fast(bboxes, nms_threshold)
        all_good_bees = filter_bees_by_NMS(all_good_bees, idxs)

    return all_good_bees

##########################################################################################

######## LOAD PARAMETERS FROM CONFIG

def load_config():
    """
        Add functionality to load previous config files
        For now, just set processing params to default ones specified in config file
    """
    import json

    processing_params = config.DEFAULT_PROCESSING_PARAMS

    return processing_params

##########################################################################################

######## SELECT VIDEO

def video_selection(dataset_dir=None):

    if dataset_dir is None:
        dataset_dir = config.DATASET_DIR

    # Select main video directory
    # ------------------------------------------------
    dataset_dirs = glob.glob(f"{dataset_dir}/*")
    print("Dataset Dirs")
    print("------------")
    for i, dataset_dir in enumerate(dataset_dirs):
        print(f"{i}. {dataset_dir}")

    user_input = int(input("Dataset dir: "))
    dataset_dir = dataset_dirs[user_input]
    print(f"'{dataset_dir}' selected.")
    # ------------------------------------------------

    # Select video from directory
    # ------------------------------------------------
    video_paths = glob.glob(f"{dataset_dir}/*")
    print("Video Paths")
    print("------------")
    for i, video_path in enumerate(video_paths):
        print(f"{i}. {video_path}")

    user_input = int(input("Video: "))
    video_path = video_paths[user_input]
    print(f"'{video_path}' selected.")

    video_name = video_path.split(os.path.sep)[-1].split('.')[0]
    return video_path, video_name

##########################################################################################

######## PREVIEW IMAGE, WITH SLICING AND CORRECTING UNEVEN ILLUMINATION

def preview_video_img(video_path, slice_r, slice_c, queen_x, queen_y, remove_background_lines=True):
    # Instantiate video streamer
    vid = utils.VideoHandler(video_path)
    for img in vid:
        img = img[slice_r, slice_c]
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        break

    if remove_background_lines:
        cleaned_img = image_processing.remove_lines(img)
        fig, ax = plt.subplots(1, 2, figsize=(14,14))
        ax[0].imshow(img, cmap=config.GRAY_CMAP)
        ax[0].set_title("Original")
        ax[1].imshow(cleaned_img, cmap=config.GRAY_CMAP)
        ax[1].set_title("Cleaned Background")
        for i in range(2):
            ax[i].scatter(queen_x, queen_y, c='r', s=120)
    else:
        fig, ax = plt.subplots(figsize=(14,14))
        ax.imshow(img, cmap=config.GRAY_CMAP)
        ax.scatter(queen_x, queen_y, c='r', s=120)

    return img, cleaned_img

##########################################################################################

######## GO THROUGH FRAMES OF VIDEO AND APPLY PROCESS_IMAGE(), OVERLAY DETECTIONS ON FRAMES

def process_video(video_path, SLICE_R, SLICE_C, processing_params, img_limit=30):
    try:
        # Instantiate video streamer
        vid = utils.VideoHandler(video_path)

        all_bees = []
        for img_i, img_BGR in enumerate(vid):
            if img_limit and img_i > img_limit:
                break

            # Slice image
            img_BGR = np.copy(img_BGR)[SLICE_R, SLICE_C]

            # Preprocess source image
            # 1. Create grayscale and, if flagged, remove lines from img
            # 2. creates overlay img source
            img_GRAY, overlay_img_src = image_processing.source_img_preprocessing(img_BGR, remove_background_lines=config.REMOVE_BACKGROUND_LINES)

            # Run image processing
            try:
                bees = process_image(img_GRAY.astype('uint8'), processing_params=processing_params)
            except Exception as e:
                print(f"\n**Error: {e}\n")
                assert False, e
            else:
                all_bees.append(bees)

            # Stdout
            sys.stdout.write(f'\rImage {img_i+1}')
            sys.stdout.flush()

    except KeyboardInterrupt:
        print("\nKeyboard Interrupt. Finished.")

    return all_bees


def process_overlay_imgs(video_path, SLICE_R, SLICE_C, QUEEN_X, QUEEN_Y,
                        all_bees, processing_params, config, img_limit=30):
    try:
        # Instantiate video streamer
        vid = utils.VideoHandler(video_path)
        overlay_imgs = []
        overlay_img_srcs = []
        for img_i, (img_BGR, bees) in enumerate(zip(vid, all_bees)):
            if img_limit and img_i > img_limit:
                break

            # Slice image
            img_BGR = np.copy(img_BGR)[SLICE_R, SLICE_C]

            # Preprocess source image
            # 1. Create grayscale and, if flagged, remove lines from img
            # 2. creates overlay img source
            img_GRAY, overlay_img_src = image_processing.source_img_preprocessing(img_BGR, remove_background_lines=config.REMOVE_BACKGROUND_LINES)
            overlay_img_srcs.append(overlay_img_src)
            overlay_img = visualization.create_overlay_img(overlay_img_src,
                                                           bees,
                                                           config=config,
                                                           processing_params=processing_params,
                                                           queen_x=QUEEN_X, queen_y=QUEEN_Y)
            overlay_imgs.append(overlay_img)

            # Stdout
            sys.stdout.write(f'\rImage {img_i+1}')
            sys.stdout.flush()

    except KeyboardInterrupt:
        print("\nKeyboard Interrupt. Finished.")

    return overlay_imgs, overlay_img_srcs

##########################################################################################

######## WRITE OVERLAYED FRAMES INTO VIDEO

def setup_video_save(video_name):
    video_dir = os.path.join(config.VID_OUTPUT_DIR, video_name)
    if not os.path.exists(video_dir):
        os.makedirs(video_dir)
        new_video_num = 0
    else:
        preexisting_vids = [ele.split(os.path.sep)[-1].replace('.mp4', '') for ele in glob.glob(f'{video_dir}/*.mp4')]
        if len(preexisting_vids) > 0:
            video_numbers = np.sort([int(ele.split("processed_")[-1]) for ele in preexisting_vids])
            max_num = max(video_numbers)
            new_video_num = max_num+1
        else:
            new_video_num = 0
    video_path = os.path.join(video_dir, f"{video_name}__processed_{new_video_num:03d}.mp4")

    return video_dir, video_path, new_video_num


def save_config(filepath, processing_params):
    import json
    with open(filepath, "w") as outfile:
        json.dump(processing_params, outfile)
