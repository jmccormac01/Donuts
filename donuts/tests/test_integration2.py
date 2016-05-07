"""Real world test of Donuts using real images from NITES
"""
import os.path
import numpy as np
from astropy.tests.helper import remote_data
from .helpers import get_test_filename
from ..donuts import Donuts


# Test the donuts code using real data
@remote_data
def test_full_integration2():
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
    test_ref_image = get_test_filename('J17490840-001.fit')
    d = Donuts(refimage=test_ref_image, image_ext=0, exposure='EXPTIME',
               normalise=True, subtract_bkg=True, prescan_width=0,
               overscan_width=0, border=64, ntiles=28)
    
    # print a summary of the setup
    d.print_summary()
    # assumes all the settings from the ref image generation
    # and calculates the shift between the images
    imlist = ['J17490840-001.fit',
              'J17490840-002.fit',
              'J17490840-003.fit']
    x_expected = [0.00, -0.73, -2.25]
    y_expected = [0.00, 2.29, 2.62]
    for image, x_ex, y_ex in zip(imlist, x_expected, y_expected):
        test_check_image = get_test_filename(image)
        result = d.measure_shift(checkimage_filename=test_check_image)
        assert np.isclose(result.x.value, x_ex, rtol=0.1, atol=0.1)
        assert np.isclose(result.y.value, y_ex, rtol=0.1, atol=0.1)
