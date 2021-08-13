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

    def __init__(self, refimage, image_ext=0, exposure='EXPTIME',
                 normalise=True, subtract_bkg=True, downweight_edges=True,
                 prescan_width=0, overscan_width=0, scan_direction='x',
                 border=64, ntiles=32, image_class=Image):
        '''Initialise and generate a reference image.
        This reference image is used for measuring frame to frame offsets.

        Parameters
        ----------
        refimage : str
            The image representing the reference frame.
        image_ext: int, optional
            The fits image extension to extract. The default is 0.
        exposure : str, optional
            Fits header keyword for exposure time. The default is `EXPTIME`.
        normalise : bool, optional
            Convert image counts to counts/s. The default is True.
        subtract_bkg : bool, optional
            Subtract the sky background. The default is True.
        downweight_edges : bool, optional
            Downweight contribution from pixels near the image edge. The default is True.
        prescan_width : int, optional
            Width of prescan region (left) in pixels. The default is 0.
        overscan_width : int, optional
            Width of overscan region (right) in pixels. The default is 0.
        scan_direction : str, optional
            Direction along which the pre/overscan regions are found ('x' | 'y')
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
        self.image_class = image_class
        self.image_ext = image_ext
        self.ntiles = ntiles
        self.exposure_keyname = exposure
        self.normalise = normalise
        self.subtract_bkg = subtract_bkg
        self.downweight_edges = downweight_edges
        self.prescan_width = prescan_width
        self.overscan_width = overscan_width
        self.scan_direction = scan_direction
        self.border = border
        self.refimage_filename = refimage

        self.reference_image = self.construct_object(self.refimage_filename)

    def construct_object(self, filename):
        '''Builds an ``image_class`` instance which performs most of the work.
        See the :class:`~donuts.image.Image` class for more information.

        Parameters
        ----------
        filename : str
            FITS file to open and build an ``image_class`` instance from.

        Returns
        -------
        ``image_class`` instance

        Raises
        ------
        None
        '''
        with fits.open(filename) as hdulist:
            hdu = hdulist[self.image_ext]
            image = hdu.data
            header = hdu.header

        image = self.image_class(image, header)
        image.preconstruct_hook()
        image.trim(
            prescan_width=self.prescan_width,
            overscan_width=self.overscan_width,
            scan_direction=self.scan_direction,
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
            if self.downweight_edges:
                image.downweight_edges()

        image.postconstruct_hook()
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
        image: :class:`~donuts.image.Image`
            Instance of an :class:`~donuts.image.Image` object, which has the ``x`` and ``y``
            atributes storing the value of the shift between the chosen image
            and reference image (passed to the :class:`~donuts.Donuts` constructor)

        Raises
        ------
        None
        '''
        checkimage = self.construct_object(checkimage_filename)
        checkimage.compute_offset(self.reference_image)
        return checkimage
