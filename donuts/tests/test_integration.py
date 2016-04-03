"""Real world test of Donuts using real images from NGTS
"""
import os.path
import numpy as np
from ..donuts import Donuts

ROOT_PATH = os.path.dirname(__file__)
DATA_DIR = os.path.join(ROOT_PATH, '..', 'data')

# Test the donuts code using real data
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
    test_ref_image = os.path.join(DATA_DIR, 'IMAGE80520160114005507.fits')
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
        test_check_image = os.path.join(DATA_DIR, image)
        x, y = d.measure_shift(checkimage=test_check_image)
        assert np.isclose(x.value, x_ex, rtol=0.1, atol=0.1)
        assert np.isclose(y.value, y_ex, rtol=0.1, atol=0.1)