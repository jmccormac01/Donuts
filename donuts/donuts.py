'''A module containing the Donuts class, used for measuring
shifts between images in CCD data.
'''
from __future__ import print_function, with_statement, division
from astropy.io import fits
from .image import Image


class Donuts(object):
    '''This class provides methods for measuring shifts between
    a series of images of the same star field. First we initialise
    the object and generate a reference image. Subsequent images are
    aligned to this frame of this reference image.

    Attributes
    ----------
    None
    '''

    def __init__(self, refimage_filename, image_ext=0, exposure_keyname='EXPTIME',
                 normalise=True, subtract_bkg=True, prescan_width=0,
                 overscan_width=0, border=64, ntiles=32):
        '''Initialise and generate a reference image.
        This reference image is used for measuring frame to frame offsets.

        Parameters
        ----------
        refimage_filename : str
            The image representing the reference frame.
        image_ext: int, optional
            The fits image extension to extract. The default is 0.
        exposure_keyname : str, optional
            Fits header keyword for exposure time. The default is `EXPTIME`.
        normalise : bool, optional
            Convert image counts to counts/s. The default is True.
        subtract_bkg : bool, optional
            Subtract the sky background. The default is True.
        prescan_width : int, optional
            Width of prescan region (left) in pixels. The default is 0.
        overscan_width : int, optional
            Width of overscan region (right) in pixels. The default is 0.
        border : int, optional
            Width of exclusion area to avoid errors from CCD edge effects.
            The default is 64.
        ntiles : int, optional
            Number of tiles used to sample the sky background.
            The default is 32.

        Returns
        -------
        None

        Raises
        ------
        None
        '''
        self.image_ext = image_ext
        self.ntiles = ntiles
        self.exposure_keyname = exposure_keyname
        self.normalise = normalise
        self.subtract_bkg = subtract_bkg
        self.prescan_width = prescan_width
        self.overscan_width = overscan_width
        self.border = border
        self.refimage_filename = refimage_filename

        self.reference_image = self.construct_object(self.refimage_filename)

    def construct_object(self, filename):
        with fits.open(filename) as hdulist:
            hdu = hdulist[self.image_ext]
            image = hdu.data
            header = hdu.header

        image = Image(image, header)
        image.trim(
            prescan_width=self.prescan_width,
            overscan_width=self.overscan_width,
            border=self.border
        )

        if self.normalise:
            image.normalise(
                exposure_keyword=self.exposure_keyname
            )

        if self.subtract_bkg:
            image.remove_background(
                ntiles=self.ntiles
            )

        image.compute_projections()

        return image

    def print_summary(self):
        '''Print a summary of the current settings

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        None
        '''
        # TODO: print more useful things!
        print('Data Summary:')
        print('\tExcluding a border of {0:d} pixels'.format(self.border))

        if self.subtract_bkg:
            print('Background Subtraction Summary:')
            print('\tUsing {0:d} x {1:d} grid of tiles'.format(self.ntiles, self.ntiles))

    def measure_shift(self, checkimage_filename):
        '''Generate a check image and measure its offset from the reference
        This is done using the same settings as the reference image.

        Parameters
        ----------
        checkimage_filename : str
            Image filename to compare to reference image

        Returns
        -------
        solution_x : float (units pixel)
            The shift required, in X, to recentre the checkimage into the reference frame
        solution_y : float (units pixel)
            The shift required, in Y, to recentre the checkimage into the reference frame

        Raises
        ------
        None
        '''
        checkimage = self.construct_object(checkimage_filename)
        checkimage.compute_offset(self.reference_image)
        return checkimage
