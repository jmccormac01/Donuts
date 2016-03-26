import numpy as np
from scipy.ndimage import gaussian_filter


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
        out[y, x] = 1.
    return out


def generate_synthetic_data(shape=(1024, 1024), background_level=100, background_sigma=5):
    '''
    Generates synthetic data, consisting of a flat gaussian background, with Gaussians for stars.
    '''
    return generate_background(
        shape=shape,
        background_level=background_level,
        background_sigma=background_sigma)
