import numpy as np
from astropy.io import fits

from donuts import Donuts

def test_no_exposure_time_keyword_uses_1(tmpdir):
    fname = tmpdir.join('test.fits')
    data = np.random.uniform(100, 1000, (1024, 1024))

    phdu = fits.PrimaryHDU(data)
    # No EXPTIME keyword
    phdu.writeto(str(fname), clobber=True)

    d = Donuts(str(fname))
    assert np.isclose(d.texp, 1.0)