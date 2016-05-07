import pytest
import numpy as np

HAS_MOCK = True
try:
    from unittest import mock
except ImportError:
    try:
        import mock
    except ImportError:
        HAS_MOCK = False

from donuts.image import Image


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


class TestImageShapeReading(object):

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

    @pytest.mark.skipif(not HAS_MOCK, reason="Cannot import the `mock` library "
            "from either the standard library, or as a package")
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
