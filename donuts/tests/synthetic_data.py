import numpy as np

def generate_background(shape, background_level, background_sigma):
    return np.random.normal(background_level, background_sigma, size=shape)

def generate_synthetic_data(shape=(1024, 1024), background_level=100, background_sigma=5):
    return generate_background(
        shape=shape,
        background_level=background_level,
        background_sigma=background_sigma)
