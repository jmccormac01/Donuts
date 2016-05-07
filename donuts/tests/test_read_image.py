import numpy as np
from astropy.tests.helper import pytest
from astropy.io import fits

from donuts import Donuts


def test_open_invalid_filename(tmpdir):
    fname = tmpdir.join('test.fits')

    with pytest.raises(IOError):
        d = Donuts(str(fname))
