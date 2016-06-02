import numpy as np
from astropy import log
from astropy import units as u
from scipy import (
    conjugate,
    polyfit,
)
from skimage.transform import resize
from scipy.fftpack import fft, ifft


class Image(object):
    '''Low level class which handles the image transformations and cross
    correlation with another ``Image`` class.
    '''

    def __init__(self, data, header=None):
        self.raw_image = data
        self.header = header or {}

        self.raw_region = None
        self.sky_background = None
        self.backsub_region = None
        self.exposure_time_value = None
        self.proj_x = None
        self.proj_y = None
        self.x = None
        self.y = None

    def normalise(self, exposure_keyword='EXPOSURE'):
        '''Convert the image data into ADU per second

        Parameters
        ----------
        exposure_keyword : str
            Fits header keyword for exposure time. The default is `EXPTIME`.

        Returns
        -------
        self : ``Image``
            The current ``Image`` instance

        Raises
        ------
        RuntimeError if the image region has not been trimmed (see the ``trim``
        method).
        '''

        try:
            self.exposure_time_value = self.header[exposure_keyword]
        except KeyError:
            log.warning(
                'Exposure time keyword "{0}" not found, assuming 1.0'.format(
                    exposure_keyword)
            )
            self.exposure_time_value = 1.0

        if self.raw_region is None:
            raise RuntimeError('Image region has not been computed.'
                               'Please ensure the `#trim` method has been called')

        self.raw_region = self.raw_region / self.exposure_time_value
        return self

    def trim(self, prescan_width=0, overscan_width=0, border=64):
        '''Remove the optional prescan and overscan from the image, as well
        as the outer `n` rows/colums of the image. Finally ensure the imaging
        region is the correct dimensions for `skimage.transform.resize`
        (i.e. a multiple of 16.)

        Parameters
        ----------
        prescan_width : int
            Remove the first ``prescan_width`` columns from the image, assuming
            the are not in the imaging region.

        overscan_width : int
            Remove the last ``overscan_width`` columns from the image, assuming
            the are not in the imaging region.

        border : int
            Ignore the first/last ``border`` rows/columns from the image,
            assuming that they are not "typical", a common case with edge
            effects in CCDs.

        Returns
        -------
        self : ``Image``
            The current ``Image`` instance

        Raises
        ------
        None


        '''
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
        '''Subtract the background from the image. See ``_generate_bkg_map`` for
        more details

        Parameters
        ----------
        ntiles : int
            Number of tiles used to sample the sky background.
            The default is 32.

        Returns
        -------
        self : ``Image``
            The current ``Image`` instance


        Raises
        ------
        None
        '''
        dim_x, dim_y = self.raw_region.shape
        tilesize_x, tilesize_y = dim_x // ntiles, dim_y // ntiles
        self.sky_background = self._generate_bkg_map(
            data=self.raw_region,
            tile_num=ntiles,
            tilesizex=tilesize_x,
            tilesizey=tilesize_y
        )
        self.backsub_region = self.raw_region - self.sky_background
        return self

    def compute_projections(self):
        '''
        Compute the projection profiles. Follows the following logic:
            * if the image has not been trimmed and not background subtracted,
            then use the raw pixels
            * otherwise if the image has been trimmed then use the trimmed
            pixels
            * otherwise use the background subtracted pixels

        See ``_projection_from_image`` for details of the projection
        calculation.

        Parameters
        ----------
        None

        Returns
        -------
        self : ``Image``
            The current ``Image`` instance

        Raises
        ------
        None
        '''
        if self.backsub_region is None and self.raw_region is None:
            region = self.raw_image
        elif self.backsub_region is None:
            region = self.raw_region
        else:
            region = self.backsub_region

        assert len(region.shape) == 2

        self.proj_x = self._projection_from_image(region, axis=0)
        self.proj_y = self._projection_from_image(region, axis=1)
        return self

    def compute_offset(self, reference_image):
        '''
        Given another ``Image`` object, compute the shift in pixel units.

        This method sets ``self.x`` and ``self.y`` to the pixel shift, and
        returns the instance so in effect the user gets a "new" ``Image``
        instance with these variables set.

        Parameters
        ----------
        reference_image : ``Image``
            The reference image to compare to. Typically when called using the
            ``Donuts`` class this will be whatever was defined as the
            "reference" image
        Returns
        -------
        self : ``Image``
            The current ``Image`` instance

        Raises
        ------
        None
        '''
        reference_image._assert_projections()
        self._assert_projections()

        results = self._cross_correlate(reference_image)
        z_pos_x, z_pos_y, phi_ref_check_m_x, phi_ref_check_m_y = results

        self.x = self._find_solution(z_pos_x, phi_ref_check_m_x)
        self.y = self._find_solution(z_pos_y, phi_ref_check_m_y)

        return self

    def preconstruct_hook(self):
        '''Hook to modify the class before any standard processing

        To add functionality, alter ``self.raw_image``
        '''
        pass

    def postconstruct_hook(self):
        '''Hook to modify the class after any standard processing

        To add functionality, alter ``self.backsub_region``
        '''
        pass

    def _assert_projections(self):
        '''Make sure the projections have been computed

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        ValueError if the projections have not been computed
        '''
        if self.proj_x is None or self.proj_y is None:
            raise ValueError(
                'Projections for %s have not been computed. '
                'Please call the #compute_projections method' % self
            )

    def _find_solution(self, z_pos, phi_ref_check_m):
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
        else:  # if z_pos[0][0] == 0:
            lra = [1, 0, -1]
            tst[0] = phi_ref_check_m[-1].real
            tst[1] = phi_ref_check_m[0].real
            tst[2] = phi_ref_check_m[1].real
            coeffs = polyfit(lra, tst, 2)
            solution = -coeffs[1] / (2 * coeffs[0])
        return solution * u.pixel

    def _cross_correlate(self, reference_image):
        # FFT of the projection spectra
        f_ref_xproj = fft(reference_image.proj_x)
        f_ref_yproj = fft(reference_image.proj_y)
        f_check_xproj = fft(self.proj_x)
        f_check_yproj = fft(self.proj_y)
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

    def _projection_from_image(self, data, axis):
        '''Function to define the actual process used to compute the
        projections.  Partially as a point to stub, perhaps as a point to
        override in a subclass, but also it is easier to understand as it's a
        pure function
        '''
        return np.sum(data, axis=axis)

    def _generate_bkg_map(self, data, tile_num, tilesizex, tilesizey):
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
        bkgmap = resize(
            coarse,
            (tilesizey * tile_num, tilesizex * tile_num),
            mode='nearest'
        )
        return bkgmap
