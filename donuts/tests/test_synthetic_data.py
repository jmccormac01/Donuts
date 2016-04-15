'''Test to run synthetic data through Donuts'''
from astropy.tests.helper import pytest
import numpy as np
from .synthetic_data import (generate_synthetic_data,
                             generate_background,
                             generate_signals,
                             save_synthetic_data,
                            )

np.random.seed(42)

def test_generate_shape():
    '''Generate test shape
    '''
    data = generate_synthetic_data(positions=[(5, 5)], shape=(10, 10))
    assert data.shape == (10, 10)

def test_background_level():
    '''Test the sky background level
    '''
    background_level = 100
    background_sigma = 0.01
    data = generate_background(
        shape=(10, 10),
        background_level=background_level,
        background_sigma=background_sigma,
    )

    expected_min = 95
    expected_max = 105
    assert expected_min < data[0][0] < expected_max


def test_generate_deltas():
    '''Test generation of star positions
    '''
    positions = [(4, 4)]
    data = generate_signals(positions=positions, shape=(10, 10))
    assert data[positions[0]] > 0.
    assert data[positions[0][1], positions[0][0] - 1] == 0.
    assert data[0, 0] == 0.


def test_generate_multiple_deltas():
    '''Test generation of multiples
    '''
    positions = [(4, 4), (8, 8), (5, 5)]
    data = generate_signals(positions=positions, shape=(10, 10))
    assert data[0, 0] == 0.
    for position in positions:
        assert data[position] > 0.
        assert data[position[1], position[0] - 1] == 0.


def test_generate_synthetic_data():
    '''Test generation of synthetic data set
    '''
    positions = [(50, 50)]
    peak_height = 1000.
    data = generate_synthetic_data(
        shape=(1024, 1024),
        background_level=100,
        background_sigma=10,
        positions=positions,
        peak_height=peak_height,
    )
    for position in positions:
        assert np.isclose(data[position], peak_height, rtol=0.8)

    assert 20. < data[0, 0] < 150


def test_write_file(tmpdir):
    '''Test outputing of fits file
    '''
    fname = tmpdir.join('out.fits')
    positions = [(50, 50)]
    peak_height = 1000.
    save_synthetic_data(
        str(fname),
        shape=(1024, 1024),
        background_level=100,
        background_sigma=10,
        positions=positions,
        peak_height=peak_height,
    )
    assert fname.isfile()


def test_no_stars_raises_error():
    '''Test that if no stars, we get an error
    '''
    positions = []
    with pytest.raises(ValueError) as err:
        data = generate_synthetic_data(positions)

    assert 'no stars' in str(err).lower()


def test_stars_outside_image_are_ignored():
    '''Test stars outside the image are ignored
    '''
    positions = [(2048, 2048)]
    signals = generate_signals((10, 10), positions)
    assert np.allclose(signals, 0.)
