import pytest
HAS_MOCK = True
try:
    from unittest import mock
except ImportError:
    try:
        import mock
    except ImportError:
        HAS_MOCK = False

from .helpers import get_test_filename
from ..donuts import Donuts

@pytest.mark.skipif(not HAS_MOCK, reason="Cannot import the `mock` library "
                    "from either the standard library, or as a package")
def test_custom_image_class():
    filename = get_test_filename('IMAGE80520160114005507.fits')
    custom_image_class = mock.MagicMock()
    d = Donuts(filename, image_class=custom_image_class)
    assert d.image_class == custom_image_class

    assert len(custom_image_class.mock_calls) > 0
