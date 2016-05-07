'''A module containing the Donuts class, used for measuring
shifts between images in CCD data.
'''
from __future__ import print_function
from __future__ import with_statement
from __future__ import division
from scipy.fftpack import ifft
from scipy.fftpack import fft
from astropy import units as u
from astropy.io import fits
from astropy import log
from scipy import conjugate
from scipy import ndimage
from scipy import polyfit
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
   
    def __cross_correlate(self, check_data):
        '''Cross correlate the reference & check images
        
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
        check_xproj = np.sum(check_data, axis=0)
        check_yproj = np.sum(check_data, axis=1)
        # FFT of the projection spectra
        f_ref_xproj = fft(self.ref_xproj)
        f_ref_yproj = fft(self.ref_yproj)
        f_check_xproj = fft(check_xproj)
        f_check_yproj = fft(check_yproj)
        # cross correlate in and look for the maximium correlation
        f_ref_xproj_conj = conjugate(f_ref_xproj)
        f_ref_yproj_conj = conjugate(f_ref_yproj)
        complex_sum_x = f_ref_xproj_conj * f_check_xproj
        complex_sum_y = f_ref_yproj_conj * f_check_yproj
        phi_ref_check_m_x = ifft(complex_sum_x)
        phi_ref_check_m_y = ifft(complex_sum_y)
        z_x = max(phi_ref_check_m_x)
        z_pos_x = np.where(phi_ref_check_m_x == z_x)
        z_y = max(phi_ref_check_m_y)
        z_pos_y = np.where(phi_ref_check_m_y == z_y)
        return z_pos_x, z_pos_y, phi_ref_check_m_x, phi_ref_check_m_y
    
    def __find_solution(self, z_pos, phi_ref_check_m):
        '''Covert the CCF into a shift solution for individual axes
        The location of the peak in the CCF is converted to a shift
        in pixels here. Sub pixel resolution is achieved by solving a
        quadratic for the minimum, using the three pixels around the peak.

        Parameters
        ----------
        z_pos : int
            The location of the peak in the CCF
        phi_ref_check_m: array-like
            The CCF array from which to extract a correction

        Returns
        -------
        solution : float
            The shift in pixels between two images along the 
            given axis

        Raises
        ------
        None
        '''
        tst = np.empty(3)
        if z_pos[0][0] <= len(phi_ref_check_m) / 2 and z_pos[0][0] != 0:
            lra = [z_pos[0][0] - 1, z_pos[0][0], z_pos[0][0] + 1]
            tst[0] = phi_ref_check_m[lra[0]].real
            tst[1] = phi_ref_check_m[lra[1]].real
            tst[2] = phi_ref_check_m[lra[2]].real
            coeffs = polyfit(lra, tst, 2)
            solution = -(-coeffs[1] / (2 * coeffs[0]))
        elif z_pos[0][0] > len(phi_ref_check_m) / 2 and z_pos[0][0] != len(phi_ref_check_m) - 1:
            lra = [z_pos[0][0] - 1, z_pos[0][0], z_pos[0][0] + 1]
            tst[0] = phi_ref_check_m[lra[0]].real
            tst[1] = phi_ref_check_m[lra[1]].real
            tst[2] = phi_ref_check_m[lra[2]].real
            coeffs = polyfit(lra, tst, 2)
            solution = len(phi_ref_check_m) + (coeffs[1] / (2 * coeffs[0]))
        elif z_pos[0][0] == len(phi_ref_check_m) - 1:
            lra = [-1, 0, 1]
            tst[0] = phi_ref_check_m[-2].real
            tst[1] = phi_ref_check_m[-1].real
            tst[2] = phi_ref_check_m[0].real
            coeffs = polyfit(lra, tst, 2)
            solution = 1 + (coeffs[1] / (2 * coeffs[0]))
        else: #if z_pos[0][0] == 0:
            lra = [1, 0, -1]
            tst[0] = phi_ref_check_m[-1].real
            tst[1] = phi_ref_check_m[0].real
            tst[2] = phi_ref_check_m[1].real
            coeffs = polyfit(lra, tst, 2)
            solution = -coeffs[1] / (2 * coeffs[0])
        return solution*u.pixel

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
        print('Data Summary:')
        print('\tIlluminated array size: {0:d} x {1:d} pixels'.format(self.dx, self.dy))
        print('\tExcluding a border of {0:d} pixels'.format(self.border))
        print('\tMeasuring shifts from central {0:d} x {1:d} pixels'.format(self.w_dimx,
                                                                            self.w_dimy))
        if self.subtract_bkg:
            print('Background Subtraction Summary:')
            print('\tUsing {0:d} x {1:d} grid of tiles'.format(self.ntiles, self.ntiles))
            print('\tEach {0:d} x {1:d} pixels'.format(self.tilesizex, self.tilesizey))

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




        # check_data = self.__get_check_data(checkimage)
        # z_pos_x, z_pos_y, phi_ref_check_m_x, phi_ref_check_m_y = self.__cross_correlate(
        #     check_data
        # )
        # solution_x = self.__find_solution(z_pos_x, phi_ref_check_m_x)
        # solution_y = self.__find_solution(z_pos_y, phi_ref_check_m_y)
        # log.debug("X: {0:.2f}".format(solution_x.value))
        # log.debug("Y: {0:.2f}".format(solution_y.value))
        # return solution_x, solution_y
