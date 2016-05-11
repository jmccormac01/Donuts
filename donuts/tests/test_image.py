import pytest
from astropy import units as u
from astropy.io import fits
from astropy.tests.helper import remote_data
import numpy as np
from .helpers import get_test_filename

try:
    from unittest import mock
except ImportError:
    import mock

from ..image import Image


def generate_image(*shape):
    return np.random.randint(1500, 2**16 - 1, size=shape).astype(
        np.uint16)


@pytest.fixture(scope='module')
def ngts_data():
    return generate_image(2048, 2048)


def compute_projections(arr):
    return np.sum(arr, axis=0), np.sum(arr, axis=1)


def test_image_stores_parameters():
    image = Image(data='a', header='b')
    assert image.raw_image == 'a'
    assert image.header == 'b'


def test_image_has_default_none_for_xy():
    image = Image(data='a', header='b')
    assert image.x is None and image.y is None


def test_header_is_optional():
    image = Image(generate_image(2048, 2048))
    assert image.header == {}


class TestTrimming(object):

    def test_without_border_or_overscans(self):
        image = Image(generate_image(2048, 2048), None)
        image.trim(border=0)
        assert image.raw_region.shape == (2048, 2048)

    def test_with_border_no_overscan(self):
        image = Image(generate_image(2048, 2048), None)
        image.trim(border=64)
        assert image.raw_region.shape == (1920, 1920)

    def test_with_border_and_overscans(self):
        image = Image(generate_image(2048, 2048), None)
        image.trim(border=64, overscan_width=20, prescan_width=20)
        assert image.raw_region.shape == (1920, 1408)

    def test_overscan_only(self):
        image = Image(generate_image(2048, 2048), None)
        image.trim(border=0, overscan_width=20)
        assert image.raw_region.shape == (2048, 1536)

    def test_prescan_only(self):
        image = Image(generate_image(2048, 2048), None)
        image.trim(border=0, prescan_width=64)
        assert image.raw_region.shape == (2048, 1984)

    def test_returns_self(self):
        image = Image(generate_image(2048, 2048))
        testvalue = image.trim()
        assert testvalue is image


class TestImageNormalisation(object):

    def test_normalise_image_with_no_exposure(self):
        image = Image(
            generate_image(2048, 2048),
            {'EXPOSURE': 1.0},
        )

        original_data = np.ones(
            (2048, 2048), dtype=np.float64) * 10.0
        image.raw_region = original_data

        image.normalise(exposure_keyword='EXPOSURE')
        assert np.allclose(image.raw_region, original_data)

    def test_normalise_image_with_exposure(self):
        image = Image(
            generate_image(2048, 2048),
            {'EXPOSURE': 10.0},
        )

        original_data = np.ones(
            (2048, 2048), dtype=np.float64) * 10.0
        image.raw_region = original_data

        image.normalise(exposure_keyword='EXPOSURE')
        assert np.allclose(image.raw_region, original_data / 10.0)

    def test_missing_exposure_time_keyword(self):
        image = Image(generate_image(2048, 2048))

        original_data = np.ones(
            (2048, 2048), dtype=np.float64) * 10.0
        image.raw_region = original_data

        image.normalise(exposure_keyword='EXPOSURE')
        assert np.allclose(image.raw_region, original_data)

    def test_no_image_region(self):
        image = Image(generate_image(2048, 2048))

        with pytest.raises(RuntimeError) as exc_info:
            image.normalise(exposure_keyword='EXPOSURE')

        assert 'image region' in str(exc_info.value).lower()

    def test_method_returns_self(self):
        image = Image(generate_image(2048, 2048))
        image.trim()
        testvalue = image.normalise()
        assert testvalue is image


class TestBackgroundSubtraction(object):

    def test_background_is_not_populated_after_init(self):
        image = Image(generate_image(2048, 2048), None)
        assert image.sky_background is None

    def test_backsub_is_not_populated_after_init(self):
        image = Image(generate_image(2048, 2048), None)
        assert image.backsub_region is None

    def test_constant_background(self):
        image = Image(generate_image(2048, 2048), None)
        image.raw_region = np.ones((2048, 2048))
        image.remove_background()
        assert np.allclose(image.sky_background, np.ones((2048, 2048)))

    # TODO: add tests for the background subtraction on non-trivial datasets

    def test_backsub_region_is_populated(self):
        image = Image(generate_image(2048, 2048), None)
        image.raw_region = np.ones((2048, 2048))
        image.remove_background()

        assert np.allclose(
            image.backsub_region,
            np.zeros((2048, 2048)))

    def test_method_returns_self(self):
        image = Image(generate_image(2048, 2048), None)
        image.raw_region = np.ones((2048, 2048))
        testvalue = image.remove_background()
        assert testvalue is image


class TestComputeProjections(object):

    def test_projections_are_not_populated_after_init(self, ngts_data):
        image = Image(data=ngts_data, header=None)
        assert image.proj_x is None
        assert image.proj_y is None

    def test_projections_are_computed_from_full_raw_image(self, ngts_data):
        image = Image(data=ngts_data, header=None)
        image.compute_projections()

        expected_proj_x, expected_proj_y = compute_projections(ngts_data)

        assert np.allclose(image.proj_x, expected_proj_x)
        assert np.allclose(image.proj_y, expected_proj_y)

    def test_projections_are_computed_from_full_raw_image(self, ngts_data):
        image = Image(data=ngts_data, header=None)
        image.compute_projections()

        assert image.proj_x.shape == (2048, )
        assert image.proj_y.shape == (2048, )

    def test_projections_with_trimmed_image(self, ngts_data):
        image = Image(data=ngts_data, header=None)
        image.trim(border=64)
        image.compute_projections()

        assert image.proj_x.shape == (1920, )
        assert image.proj_y.shape == (1920, )

    def test_projections_with_backsub_image(self):
        stub_data = np.ones((2048, 2048))
        image = Image(data=stub_data, header=None)
        image.trim(border=64)
        image.remove_background()

        background_image = image.backsub_region

        with mock.patch.object(image, '_projection_from_image') as proj_fn:
            image.compute_projections()

        calls = [
            mock.call(background_image, axis=0),
            mock.call(background_image, axis=1),
        ]
        proj_fn.assert_has_calls(calls)

    def test_method_returns_self(self):
        image = Image(generate_image(2048, 2048))
        image.trim()
        testvalue = image.compute_projections()
        assert testvalue is image


class TestComputeOffset(object):

    @staticmethod
    def build_image(image_filename):
        with fits.open(image_filename) as infile:
            image = Image(infile[0].data, infile[0].header)

        image.trim()
        image.normalise()
        image.remove_background()
        image.compute_projections()
        return image

    @staticmethod
    def is_close(x, y):
        if isinstance(x, u.Quantity):
            x = x.value
        if isinstance(y, u.Quantity):
            y = y.value

        return np.isclose(x, y, rtol=0.1, atol=0.1)

    @remote_data
    def test_same_object_gives_same_offset(self, ngts_data):
        i1 = self.build_image(
            get_test_filename('IMAGE80520160114005507.fits')
        )
        i2 = self.build_image(
            get_test_filename('IMAGE80520160114005507.fits')
        )
        i2.compute_offset(i1)

        assert self.is_close(i2.x, 0.)
        assert self.is_close(i2.y, 0.)

    def test_error_without_projection_calculation(self, ngts_data):
        image = Image(ngts_data)
        assert image.proj_x is None and image.proj_y is None

        with pytest.raises(ValueError) as exc_info:
            image.compute_offset(image)

        assert 'Please call the #compute_projections method' in str(exc_info.value)

    @remote_data
    def test_known_offsets(self):
        ref_image = self.build_image(
            get_test_filename('IMAGE80520160114005507.fits')
        )
        test_image = self.build_image(
            get_test_filename('IMAGE80520160114005520.fits')
        )

        expected = (-0.09, 0.24)

        test_image.compute_offset(ref_image)

        assert self.is_close(test_image.x, expected[0])
        assert self.is_close(test_image.y, expected[1])

    # Lots of mocking for functions that will fail without proper setup
    @mock.patch.object(Image, '_assert_projections')
    @mock.patch.object(Image, '_cross_correlate',
                       return_value=(None, None, None, None))
    @mock.patch.object(Image, '_find_solution')
    def test_method_returns_self(self,
                                 assert_projections,
                                 cross_correlate,
                                 find_solution,
                                 ):
        image = Image(generate_image(2048, 2048))
        image2 = Image(generate_image(2048, 2048))

        testvalue = image2.compute_offset(image)
        assert testvalue is image2

def test_no_exposure_time_keyword_uses_1():
    data = np.random.uniform(100, 1000, (1024, 1024))

    image = Image(data)
    image.trim()
    image.normalise()

    assert np.isclose(image.exposure_time_value, 1.0)

def test_image_class_has_pre_setup_hook():
    image = Image(generate_image(2048, 2048))
    assert hasattr(image, 'preconstruct_hook')

def test_image_class_has_post_setup_hook():
    assert hasattr(Image, 'postconstruct_hook')
