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

# to do:
#   add try/except for opening images etc - IN PROGRESS
#   set flag when reference image is made, do not run check if no reference - IN PROGRESS

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
                 normalise=True, subtract_bkg=True, prescan_width=0,
                 overscan_width=0, border=64, ntiles=32):
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
        self.normalise = normalise
        self.refimage = refimage
        self.subtract_bkg = subtract_bkg
        self.prescan_width = prescan_width
        self.overscan_width = overscan_width
        self.border = border
        self.solution_x = 0.0
        self.solution_y = 0.0
        self.base = 512

        # define here first or pylint cries
        self.checkimage = None
        self.check_data = None
        self.check_image_section = None
        self.check_bkgmap = None
        self.check_xproj = None
        self.check_yproj = None

        with fits.open(refimage) as h:
            try:
                exposure_time_value = h[self.image_ext].header[exposure]
            except KeyError:
                log.warning('Exposure time keyword "{0}" not found, assuming 1.0'.format(
                    exposure))
                exposure_time_value = 1.0

            self.texp = float(exposure_time_value)
            # get image dimmensions
            if self.overscan_width != 0:
                self.image_section = h[self.image_ext].data[:, self.prescan_width:-self.overscan_width]
            else:
                self.image_section = h[self.image_ext].data[:, self.prescan_width:]
            self.dy, self.dx = self.image_section.shape
            self.isRefImage =  True
        
        # check if the CCD is a funny shape. Normal CCDs should divide by 16 with
        # no remainder. NITES for example does not (1030,1057) instead of (1024,1024)
        rx = self.dx % 16
        ry = self.dy % 16
        if rx > 0 or ry > 0:
            self.dimx = int(self.dx // self.base) * self.base
            self.dimy = int(self.dy // self.base) * self.base
        else:
            self.dimx = self.dx
            self.dimy = self.dy
        # get the reference data, with tweaked shape if needed
        self.ref_data = self.image_section[self.border:self.dimy - self.border,
                                           self.border:self.dimx - self.border]
        # get the working image dimensions after removing the border
        self.w_dimy, self.w_dimx = self.ref_data.shape
        # set up tiles for bkg subtract
        self.tilesizex = self.w_dimx // self.ntiles
        self.tilesizey = self.w_dimy // self.ntiles
        # adjust image if requested
        if self.subtract_bkg:
            self.bkgmap = self.__generate_bkg_map(self.ref_data, self.ntiles,
                                                  self.tilesizex, self.tilesizey)
            self.ref_data = self.ref_data - self.bkgmap
        if self.normalise:
            self.ref_data = self.ref_data / self.texp

        self.ref_xproj = np.sum(self.ref_data, axis=0)
        self.ref_yproj = np.sum(self.ref_data, axis=1)

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

    def __get_check_data(self, checkimage):
        '''Grab the check_data array
        This is done using the same settings as the reference image

        Parameters
        ----------
        checkimage : str
            Image filename to compare to reference image

        Returns
        -------
        None

        Raises
        ------
        None
        '''
        self.checkimage = checkimage
        with fits.open(self.checkimage) as h:
            if self.overscan_width != 0:
                self.check_image_section = h[self.image_ext].data[:, self.prescan_width:-self.overscan_width]
            else:
                self.check_image_section = h[self.image_ext].data[:, self.prescan_width:]
            self.check_data = self.check_image_section[self.border:self.dimy - self.border,
                                                       self.border:self.dimx - self.border]
        # adjust image if requested - same as reference
        if self.subtract_bkg:
            self.check_bkgmap = self.__generate_bkg_map(self.check_data, self.ntiles,
                                                        self.tilesizex, self.tilesizey)
            self.check_data = self.check_data - self.check_bkgmap
        if self.normalise:
            self.check_data = self.check_data / self.texp

    def __cross_correlate(self):
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
        self.check_xproj = np.sum(self.check_data, axis=0)
        self.check_yproj = np.sum(self.check_data, axis=1)
        # FFT of the projection spectra
        f_ref_xproj = fft(self.ref_xproj)
        f_ref_yproj = fft(self.ref_yproj)
        f_check_xproj = fft(self.check_xproj)
        f_check_yproj = fft(self.check_yproj)
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
        return solution

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

    def measure_shift(self, checkimage):
        '''Generate a check image and measure its offset from the reference
        This is done using the same settings as the reference image.

        Parameters
        ----------
        checkimage : str
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
        self.__get_check_data(checkimage)
        z_pos_x, z_pos_y, phi_ref_check_m_x, phi_ref_check_m_y = self.__cross_correlate()
        self.solution_x = self.__find_solution(z_pos_x, phi_ref_check_m_x)
        self.solution_y = self.__find_solution(z_pos_y, phi_ref_check_m_y)
        log.debug("X: {0:.2f}".format(self.solution_x))
        log.debug("Y: {0:.2f}".format(self.solution_y))
        return self.solution_x*u.pixel, self.solution_y*u.pixel
