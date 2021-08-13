'''Synthetic data test suite'''
import numpy as np
import astropy
from astropy.io import fits
from scipy.ndimage import gaussian_filter

astropy_major_version = int(astropy.__version__.split('.')[0])

def generate_background(shape, background_level, background_sigma):
    '''
    Generates a Gaussian background signal, with the supplied level and sigma
    '''
    return np.random.normal(background_level, background_sigma, size=shape)


def generate_signals(shape, positions):
    '''
    Generates delta functions at the specified positions
    '''
    out = np.zeros(shape)
    for position in positions:
        x, y = position
        try:
            out[y, x] = 1.
        except IndexError:
            pass
    return out


def generate_synthetic_data(
        positions,
        shape=(1024, 1024),
        background_level=100,
        background_sigma=5,
        peak_height=1000.,
        psf_pix=1.6,
):
    '''
    Generates synthetic data, consisting of a flat gaussian background, with Gaussians for stars.
    '''
    background = generate_background(
        shape=shape,
        background_level=background_level,
        background_sigma=background_sigma)
    deltas = generate_signals(shape=shape, positions=positions)
    if np.allclose(deltas, 0.):
        raise ValueError("No stars found in the synthetic image")
    blurred = gaussian_filter(deltas, sigma=2.35 * psf_pix)
    blurred_max = blurred.max()
    blurred_scaled = blurred * peak_height / blurred_max
    return background + blurred_scaled


def save_synthetic_data(filename, *args, **kwargs):
    '''Save the synthetic data to file '''
    data = generate_synthetic_data(*args, **kwargs)
    phdu = fits.PrimaryHDU(data)
    if astropy_major_version < 2:
        phdu.writeto(filename, clobber=True)
    else:
        phdu.writeto(filename, overwrite=True)
