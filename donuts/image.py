import numpy as np


class Image(object):
    '''Encapsulate the transformations applied to images
    '''

    def __init__(self, data, header):
        self.raw_image = data
        self.raw_region = None
        self.sky_background = None
        self.backsub_region = None
        self.header = header
        self.proj_x = None
        self.proj_y = None
        self.x = None
        self.y = None

    def trim(self, prescan_width=0, overscan_width=0, border=64):
        if overscan_width > 0 and prescan_width > 0:
            image_section = self.raw_image[:, prescan_width:-overscan_width]
        elif overscan_width > 0:
            image_section = self.raw_image[:, :-overscan_width]
        elif prescan_width > 0:
            image_section = self.raw_image[:, prescan_width:]
        else:
            image_section = self.raw_image

        dy, dx = image_section.shape  

        # check if the CCD is a funny shape. Normal CCDs should divide by 16
        # with no remainder. NITES for example does not (1030,1057)
        # instead of (1024,1024)
        rx = dx % 16
        ry = dy % 16
        base = 512

        if rx > 0 or ry > 0:
            dimx = int(dx // base) * base
            dimy = int(dy // base) * base
        else:
            dimx = dx
            dimy = dy

        # get the reference data, with tweaked shape if needed
        self.raw_region = image_section[
            border:dimy - border, border:dimx - border
        ]
