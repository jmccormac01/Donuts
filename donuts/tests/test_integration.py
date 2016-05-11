"""Real world test of Donuts using real images from NGTS
"""
import os.path
import pytest
import numpy as np
from astropy.tests.helper import remote_data
from .helpers import get_test_filename
from ..donuts import Donuts
from ..image import Image


# Test the donuts code using real data
@remote_data
def test_full_integration():
    """Test real world usage of Donuts with real data

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # initialise the class with the settings needed
    # upon generation of the reference image

    test_ref_image = get_test_filename('IMAGE80520160114005507.fits')
    d = Donuts(refimage=test_ref_image, image_ext=0, exposure='EXPOSURE',
               normalise=True, subtract_bkg=True, prescan_width=20,
               overscan_width=20, border=64, ntiles=32)
    # print a summary of the setup
    d.print_summary()
    # assumes all the settings from the ref image generation
    # and calculates the shift between the images
    imlist = ['IMAGE80520160114005507.fits',
              'IMAGE80520160114005520.fits',
              'IMAGE80520160114005533.fits']
    x_expected = [0.00, -0.09, 0.01]
    y_expected = [0.00, 0.24, 0.14]
    for image, x_ex, y_ex in zip(imlist, x_expected, y_expected):
        test_check_image = get_test_filename(image)

        result = d.measure_shift(checkimage_filename=test_check_image)
        assert np.isclose(result.x.value, x_ex, rtol=0.1, atol=0.1)
        assert np.isclose(result.y.value, y_ex, rtol=0.1, atol=0.1)


@remote_data
def test_custom_preconstruct_code():
    test_ref_image = get_test_filename('IMAGE80520160114005507.fits')

    class CustomImage(Image):

        def preconstruct_hook(self):
            raise TypeError('bad')

    with pytest.raises(TypeError):
        d = Donuts(refimage=test_ref_image, image_ext=0,
                   exposure='EXPOSURE', normalise=True,
                   subtract_bkg=True, prescan_width=20, overscan_width=20,
                   border=64, ntiles=32, image_class=CustomImage)


@remote_data
def test_custom_postconstruct_code():
    test_ref_image = get_test_filename('IMAGE80520160114005507.fits')

    class CustomImage(Image):

        def postconstruct_hook(self):
            raise TypeError('bad')

    with pytest.raises(TypeError):
        d = Donuts(refimage=test_ref_image, image_ext=0,
                   exposure='EXPOSURE', normalise=True,
                   subtract_bkg=True, prescan_width=20, overscan_width=20,
                   border=64, ntiles=32, image_class=CustomImage)
