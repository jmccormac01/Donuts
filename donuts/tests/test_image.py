import pytest
import numpy as np

from donuts.image import Image


def generate_image(shape):
    return np.random.randint(1500, 2**16 - 1, size=shape).astype(
        np.uint16)


@pytest.fixture(scope='module')
def ngts_data():
    return generate_image((2048, 2048))


def compute_projections(arr):
    return np.sum(arr, axis=0), np.sum(arr, axis=1)


def test_image_stores_parameters():
    image = Image(data='a', header='b')
    assert image.raw_image == 'a'
    assert image.header == 'b'


def test_image_has_default_none_for_xy():
    image = Image(data='a', header='b')
    assert image.x is None and image.y is None


class TestImageShapeReading(object):

    def test_without_border_or_overscans(self):
        image = Image(generate_image((2048, 2048)), None)
        image.trim(border=0)
        assert image.raw_region.shape == (2048, 2048)

    def test_with_border_no_overscan(self):
        image = Image(generate_image((2048, 2048)), None)
        image.trim(border=64)
        assert image.raw_region.shape == (1920, 1920)

    def test_with_border_and_overscans(self):
        image = Image(generate_image((2048, 2048)), None)
        image.trim(border=64, overscan_width=20, prescan_width=20)
        assert image.raw_region.shape == (1920, 1408)

    def test_overscan_only(self):
        image = Image(generate_image((2048, 2048)), None)
        image.trim(border=0, overscan_width=20)
        assert image.raw_region.shape == (2048, 1536)

    def test_prescan_only(self):
        image = Image(generate_image((2048, 2048)), None)
        image.trim(border=0, prescan_width=64)
        assert image.raw_region.shape == (2048, 1984)
