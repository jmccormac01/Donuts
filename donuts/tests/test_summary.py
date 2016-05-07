import pytest
from astropy.io import fits
from .helpers import get_test_filename
from ..donuts import Donuts


def build_object(subtract_background):
    filename = get_test_filename('IMAGE80520160114005507.fits')
    return Donuts(
        refimage=filename,
        exposure='EXPOSURE',
        prescan_width=20,
        overscan_width=20,
        ntiles=32,
        border=64,
        subtract_bkg=subtract_background
    )


@pytest.mark.parametrize('expected', [
    'Data Summary',
    'Excluding a border of 64 pixels',
    'Using 32 x 32 grid of tiles',
])
def test_printing(expected, capsys):
    dobj = build_object(subtract_background=True)
    dobj.print_summary()
    out, _ = capsys.readouterr()
    assert expected in out


def test_with_background_subtraction(capsys):
    dobj = build_object(subtract_background=True)
    dobj.print_summary()
    out, _ = capsys.readouterr()
    assert 'Background Subtraction Summary' in out


def test_without_background_subtraction(capsys):
    dobj = build_object(subtract_background=False)
    dobj.print_summary()
    out, _ = capsys.readouterr()
    assert 'Background Subtraction Summary' not in out
