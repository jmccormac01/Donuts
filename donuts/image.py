import numpy as np
import warnings
from astropy import log
from scipy import ndimage


class Image(object):
    '''Encapsulate the transformations applied to images
    '''

    def __init__(self, data, header=None):
        self.raw_image = data
        self.header = header or {}

        self.raw_region = None
        self.sky_background = None
        self.backsub_region = None
        self.proj_x = None
        self.proj_y = None
        self.x = None
        self.y = None

    def normalise(self, exposure_keyword='EXPOSURE'):
        try:
            exposure_time_value = self.header[exposure_keyword]
        except KeyError:
            log.warning(
                'Exposure time keyword "{0}" not found, assuming 1.0'.format(
                    exposure_keyword)
            )
            exposure_time_value = 1.0

        if self.raw_region is None:
            raise RuntimeError('Image region has not been computed.'
                               'Please ensure the `#trim` method has been called')

        self.raw_region = self.raw_region / exposure_time_value
        return self

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
        return self

    def remove_background(self, ntiles=32):
        dim_x, dim_y = self.raw_region.shape
        tilesize_x, tilesize_y = dim_x // ntiles, dim_y // ntiles
        self.sky_background = self.__generate_bkg_map(
            data=self.raw_region,
            tile_num=ntiles,
            tilesizex=tilesize_x,
            tilesizey=tilesize_y
        )

    def __generate_bkg_map(self, data, tile_num, tilesizex, tilesizey):
        '''Create a background map.
        This map may be subtracted from each image before doing the cross correlation

        Parameters
        ----------
        data : array-like
            Image array from which to measure sky background
        tile_num : int
            Number of tiles along each axis
        tilesizex : int
            Size of tiles in X, pixels
        tilesizey : int
            Size of tiles in Y, pixels

        Returns
        -------
        bkgmap : array-like
            2D map of sky background. This can then be subtracted
            from each image to improve the shift measurement

        Raises
        ------
        None
        '''
        # create coarse map
        coarse = np.empty((tile_num, tile_num))
        for i in range(0, tile_num):
            for j in range(0, tile_num):
                coarse[i][j] = np.median(data[(i * tilesizey):(i + 1) * tilesizey,
                                              (j * tilesizex):(j + 1) * tilesizex])
        # resample it out to data size
        bkgmap = ndimage.zoom(coarse, (tilesizey, tilesizex), order=2)
        return bkgmap
