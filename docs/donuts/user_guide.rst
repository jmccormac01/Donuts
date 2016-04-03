**********
User guide
**********

Donuts computes the translational shift between two images of the same
star field. It has been used with great success for multiple instruments
including defocused observations.

In general the user must supply one image as the "reference image",
which is treated as 0,0 in the translational coordinate system. All
other images are then aligned to this image.

Example usage
-------------

The example below shows typical usage of the ``donuts`` api::

    >>> from donuts import Donuts
    >>> reference_image_name = 'refimage.fits'
    >>> science_image_names = ['image1.fits', 'image2.fits']
    >>> # Construct a donuts object
    >>> d = Donuts(
    ...   refimage=reference_image_name,
    ...   exposure='EXPOSURE',
    ...   overscan_width=20,
    ...   prescan_width=20)
    >>> # for each image, compute the x/y translation required
    ... # to align the images onto the reference image
    >>> for image in science_image_names:
    ...     x, y = d.measure_shift(checkimage=image)
    ...     print(x, y)

Further detail
--------------

More explanation of the ``Donuts`` constructor parameters are given
below.

Defining the reference image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following section describes some settings that can be applied to
each image. It is important that the process applied to the reference
image is applied to each science image to ensure that only the field
position changes between images.

Image extension
```````````````

The image extension for each image is set with the ``image_ext``
parameter, which defaults to 0 (the primary HDU). If your image is in
another extension, specify it here.

Normalisation
`````````````

Image normalisation (by the exposure time) can be enabled through the
``normalise`` parameter, and the exposure time chosen with the
``exposure`` parameter, which defaults to ``EXPTIME``. By default this
feature is on.

Background subtraction
``````````````````````

It is often useful to subtract the sky background from each image,
leaving just the stars to match against. The sky background is computed
on a coarse grid across the image and then interpolated between to
compute the per-pixel sky background level.

Prescan and overscan removal
````````````````````````````

If your images contain either prescan or overscan regions, set the width
of each using the ``prescan_width`` and ``overscan_width`` values (both
default to 0).

Border
``````

To remove the edge N pixels of the image before aligning (often useful),
use the ``border`` parameter. This defaults to 64.

Computing per-image offsets
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With a ``Donuts`` object, the x/y translation required to align the
image into the reference image is given by the ``measure_shift`` method,
e.g.::

    >>> x, y = d.measure_shift('image.fits')

x and y are given in pixel units.
