'''A module containing the Donuts class, used for measuring
shifts between images in CCD data.
'''
from __future__ import print_function, with_statement, division
from astropy.io import fits
import numpy as np
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
                 border=64, ntiles=32, calculation_area_override=None,
                 image_pixel_mask=None, image_class=Image):
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
        calculation_area_override : list|tuple, optional
            Manually supplied coordinates for shift calculation image area
            Image region is defined as (lower_y, upper_y, lower_x, upper_x)
            e.g. to calculate shifts using the lower left corner with a 500 pix
            square we'd supply:
                (0, 500, 0, 500)
        image_pixel_mask : array | str, optional
            Array of booleans (0|1 or False|True) where the affirmative corresponds
            to the locations of pixels to be masked out from shift calculations
            (e.g. the location of hot-pixels). This boolean array must have the same
            shape as the imager sensor array, including any pre/overscan areas. The mask
            is applied to the untrimmed image immediately after loading the data.

            If a str is supplied this is assumed to be the path to a fits image on
            disc that contains the image mask in the first (0th) image extension.
            As above, the fits image must be a boolean array (0|1 or False|True) of
            the same shape at the imager, including any pre/overscan regions.

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

        # store the image geometry region corners
        # NOTE: ugh, adding manual image area selection adds a crazy number of checks
        # to be made, e.g. that upper bounds are greater than lower bounds, that they
        # are positive, that they fall within the image etc etc.
        # This is an advanced feature and people using this should know enough to
        # supply the values correctly. I'll add a bunch of error checking here if
        # it becomes an issue....
        if calculation_area_override:
            self.image_cly = calculation_area_override[0]
            self.image_cuy = calculation_area_override[1]
            self.image_clx = calculation_area_override[2]
            self.image_cux = calculation_area_override[3]
            self.image_geometry_set = True
        else:
            self.image_cly = None
            self.image_cuy = None
            self.image_clx = None
            self.image_cux = None
            self.image_geometry_set = False

        # determine if we're using a mask, if so is it a str (load a fits file) or
        # an array which we can directly apply
        if image_pixel_mask is not None and isinstance(image_pixel_mask, str):
            # load the fits image
            with fits.open(image_pixel_mask) as fitsfile:
                self.image_pixel_mask = fitsfile[0].data
        elif image_pixel_mask is not None and isinstance(image_pixel_mask, np.ndarray):
            # load the array as the mask
            self.image_pixel_mask = image_pixel_mask
        elif image_pixel_mask is None:
            # no mask, leave it as None
            self.image_pixel_mask = image_pixel_mask
        else:
            type_err_str = """Unhandled mask type, please supply one of the following:
              1: A numpy array of booleans
              2: A str containing the path to a fits image with the mask
              3: None, for no masking"""
            raise TypeError(type_err_str)

        # make a reference image object
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
        # load the data from disc
        with fits.open(filename) as hdulist:
            hdu = hdulist[self.image_ext]
            image = hdu.data
            header = hdu.header

        # use masked arrays regardless of masking or not
        # by this point we should have an array or None
        # None will work, but an array of the wrong shape will throw an error
        # print some info if the shapes mismatch
        try:
            masked_image = np.ma.array(image, mask=self.image_pixel_mask, fill_value=0)
        except np.ma.core.MaskError:
            mask_err_str = """Wrong mask shape for image
            Image shape: {}
            Mask shape: {}""".format(image.shape, self.image_pixel_mask.shape)
            raise Exception(mask_err_str)

        # create the image object
        image = self.image_class(masked_image, header)

        # run the preconstruct hook
        image.preconstruct_hook()

        # get the image geometry
        if not self.image_geometry_set:
            cly, cuy, clx, cux = image.calculate_image_geometry(
                prescan_width=self.prescan_width,
                overscan_width=self.overscan_width,
                scan_direction=self.scan_direction,
                border=self.border,
                ntiles=self.ntiles
            )
            self.image_cly = cly
            self.image_cuy = cuy
            self.image_clx = clx
            self.image_cux = cux
            self.image_geometry_set = True

        # trim the image to match the image geometry
        image.trim(self.image_cly,
                   self.image_cuy,
                   self.image_clx,
                   self.image_cux)

        # apply exposure time normalisation
        if self.normalise:
            image.normalise(
                exposure_keyword=self.exposure_keyname
            )

        # remove the sky background
        if self.subtract_bkg:
            image.remove_background(
                ntiles=self.ntiles
            )
            if self.downweight_edges:
                image.downweight_edges()

        # run the postconstruct hook
        image.postconstruct_hook()

        # compute the x and y image projections
        image.compute_projections()

        # return the processed image object, ready for shift calculations
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
