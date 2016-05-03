from donuts.image import Image


def test_image_stores_parameters():
    image = Image(data='a', header='b')
    assert image.raw_image == 'a'
    assert image.header == 'b'
