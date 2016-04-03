**********
User guide
**********

Donuts computes the translational shift between two images of the same
star field. It has been used with great success for multiple instruments
including defocused observations.

In general the user must supply one image as the "reference image",
which is treated as 0,0 in the translational coordinate system. All
other images are then aligned to this image.

Defining the reference image
----------------------------

The following section describes some settings that can be applied to
each image. It is important that the process applied to the reference
image is applied to each science image to ensure that only the field
position changes between images.

Image extension
~~~~~~~~~~~~~~~

The image extension for each image is set with the ``image_ext``
parameter, which defaults to 0 (the primary HDU). If your image is in
another extension, specify it here.

Normalisation
~~~~~~~~~~~~~

Image normalisation (by the exposure time) can be enabled through the
``normalise`` parameter, and the exposure time chosen with the
``exposure`` parameter, which defaults to ``EXPTIME``. By default this
feature is on.

Background subtraction
~~~~~~~~~~~~~~~~~~~~~~

It is often useful to subtract the sky background from each image,
leaving just the stars to match against. The sky background is computed
on a coarse grid across the image and then interpolated between to
compute the per-pixel sky background level.

Prescan and overscan removal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

e

Border
~~~~~~

f



Computing per-image offsets
---------------------------
