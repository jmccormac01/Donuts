import numpy as np


def generate_synthetic_data(shape=(1024, 1024), background_level=100, background_sigma=5):
    background = np.random.normal(background_level, background_sigma, size=shape)
    return background
