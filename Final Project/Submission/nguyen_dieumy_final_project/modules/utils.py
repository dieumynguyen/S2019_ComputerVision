import cv2
import matplotlib.pyplot as plt

class VideoHandler:
    def __init__(self, vid_path, color=True, img_limit=None, img_skip=1, start_i=0, end_i=None):
        self.vid_path = vid_path
        self.img_limit = img_limit
        self.img_skip = img_skip
        self.start_i = start_i
        self.end_i = end_i
        
        self._open_stream(vid_path)
        self.num_images_loaded = 0
        
    def _open_stream(self, vid_path):
        self.cap = cv2.VideoCapture(vid_path)
        
    def __iter__(self):
        self.frame_i = 0
        self.num_images_loaded = 0
        return self

    def __next__(self):
        # Read frame and increment frame counter
        ret, frame = self.cap.read()
        self.frame_i += 1
        
        # Check for image limit
        condition_1 = self.img_limit and self.num_images_loaded >= self.img_limit
        condition_2 = self.end_i is not None and self.frame_i >= self.end_i
        if condition_1 or condition_2:
            raise StopIteration
        # Check image skip
        elif (self.frame_i % self.img_skip != 0) or (self.frame_i < self.start_i):
            frame = self.__next__()
        else:
            if frame is None:
                raise StopIteration
            
            self.num_images_loaded += 1
                
        return frame
    
def imgs2vid(imgs, outpath, frame_interval=1, fps=12):
    imgs = imgs[::frame_interval]
    height, width = imgs[0].shape[0:2]
        
    fourcc = cv2.VideoWriter_fourcc("m", "p", "4", "v")
    video = cv2.VideoWriter(outpath, fourcc, fps, (width, height), True)
    
    for img in imgs:
        video.write(img)
    
    cv2.destroyAllWindows()
    video.release()
 