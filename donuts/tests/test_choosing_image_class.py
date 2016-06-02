import pytest
try:
    from unittest import mock
except ImportError:
    import mock

from astropy.tests.helper import remote_data

from .helpers import get_test_filename
from ..donuts import Donuts

@remote_data
def test_custom_image_class():
    filename = get_test_filename('IMAGE80520160114005507.fits')
    custom_image_class = mock.MagicMock()
    d = Donuts(filename, image_class=custom_image_class)
    assert d.image_class == custom_image_class

    assert len(custom_image_class.mock_calls) > 0
