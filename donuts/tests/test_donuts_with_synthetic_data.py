import pytest
import numpy as np
from astropy.io import fits
from .synthetic_data import (
    generate_synthetic_data,
    save_synthetic_data,
)
from ..donuts import Donuts

np.random.seed(42)

OVERSCAN_WIDTH = 20
PRESCAN_WIDTH = 20


def write_data(filename, data):
    full_array = np.zeros((data.shape[0], data.shape[1] + OVERSCAN_WIDTH + PRESCAN_WIDTH))
    full_array[:, PRESCAN_WIDTH:-OVERSCAN_WIDTH] = data

    phdu = fits.PrimaryHDU(full_array)
    phdu.header['EXPOSURE'] = 1.
    phdu.writeto(filename, clobber=True)


def shift_positions(positions, dx=0, dy=0):
    out = []
    for x, y in positions:
        out.append((x + dx, y + dy))
    return out


def test_same_image_gives_0_offsets(tmpdir):
    refimage = tmpdir.join('refimage.fits')
    scienceimage = tmpdir.join('scienceimage.fits')

    nstars = 500
    positions = [row for row in zip(
        np.random.uniform(0, 1023, nstars),
        np.random.uniform(0, 1023, nstars))
    ]

    data = generate_synthetic_data(positions)
    write_data(str(refimage), data)
    write_data(str(scienceimage), data)

    d = Donuts(refimage=str(refimage), image_ext=0, exposure='EXPOSURE',
               normalise=True, subtract_bkg=True, prescan_width=PRESCAN_WIDTH,
               overscan_width=OVERSCAN_WIDTH, boarder=8, ntiles=16)
    x, y = d.measure_shift(str(scienceimage))
    assert np.isclose(x, 0.) and np.isclose(y, 0.)


@pytest.mark.parametrize('dx,dy', [
    (0., 1.),
    (-1., 1.),
    (-5., -0.2),
    (2., 2.),
])
def test_known_offset(dx, dy, tmpdir):
    refimage = tmpdir.join('refimage.fits')
    scienceimage = tmpdir.join('scienceimage.fits')

    nstars = 1500
    refimage_positions = [row for row in zip(
        np.random.uniform(0, 1023, nstars),
        np.random.uniform(0, 1023, nstars))
    ]
    science_image_positions = shift_positions(
        refimage_positions, dx=dx, dy=dy
    )

    refimage_data = generate_synthetic_data(refimage_positions)
    scienceimage_data = generate_synthetic_data(science_image_positions)
    write_data(str(refimage), refimage_data)
    write_data(str(scienceimage), scienceimage_data)

    d = Donuts(refimage=str(refimage), image_ext=0, exposure='EXPOSURE',
               normalise=True, subtract_bkg=True, prescan_width=PRESCAN_WIDTH,
               overscan_width=OVERSCAN_WIDTH, boarder=8, ntiles=16)
    x, y = d.measure_shift(str(scienceimage))

    # Fairly large margin for the tolerence as the values are relatively uncertain
    assert np.isclose(x, -dx, rtol=0.5, atol=0.1)
    assert np.isclose(y, -dy, rtol=0.5, atol=0.1)
