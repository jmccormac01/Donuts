**********
User guide
**********

Donuts computes the translational shift between two images of the same
star field. It has been used with great success for multiple instruments
including defocused observations.

To measure the shifts present in a series of images requires the
definition of a reference frame. An image must be chosen from the series
to be the reference. This is treated as 0,0 in the translational 
coordinate syst. Donuts then measures the shift of each subsequent
image with respect to the reference. Reference images should be free 
from unwanted artefacts (satellite trails etc.)

Example usage
-------------

The example below shows typical usage of the ``donuts`` api:

.. doctest-skip::

    >>> from donuts import Donuts
    >>> reference_image_name = 'image1.fits'
    >>> science_image_names = ['image2.fits', 'image2.fits']
    >>> # Construct a donuts object
    >>> d = Donuts(
    ...   refimage=reference_image_name,
    ...   image_ext=0,
    ...   overscan_width=20,
    ...   prescan_width=20,
    ...   border=64,
    ...   normalise=True,
    ...   exposure='EXPOSURE',
    ...   subtract_bkg=True,
    ...   ntiles=32)
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
each image. It is important that the process applied to 
the reference image is applied to each science image to ensure that 
only the field position changes between images, hence the settings 
chosen during initialization are used for each subseqient image.

Image extension
```````````````

The fits image extension to be used in the calculation is set with 
the ``image_ext`` parameter, which defaults to 0 (the primary HDU). 
If your image is in another extension, specify it here. 


Prescan and overscan removal
````````````````````````````

If your images contain either prescan or overscan regions, set the width
of each using the ``prescan_width`` and ``overscan_width`` values (both
default to 0).


Recommended image corrections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each correction below is optional but is recommended. If you find that 
an extra correction is required please submit an issue on our 
`Github Issues <https://github.com/jmccormac01/Donuts/issues>`_ page.


Border
``````

Border pixels typically have a non-linear response and are recommended to be
excluded. To remove the edge N pixels of the image before aligning,
use the ``border`` parameter. This defaults to 64. Set this parameter to 
0 if you wish to ignore the border exclusion. 

Exposure time normalisation
```````````````````````````

Exposure time normalisation can be enabled through the
``normalise`` parameter, and the exposure time chosen with the
``exposure`` parameter, which defaults to ``EXPTIME``. By default this
feature is on. This allows for varying exposure times without the need for 
separate analysis. 

Background subtraction
``````````````````````

By design, Donuts will lock onto the most dominant 'feature' in a series
of images. Under normal circumstances this will be the forrest of stellar 
profiles in the collapsed image projections. However, under new moon 
conditions for example, a strong vignetting pattern might be the dominant 
feature and will result in poor performance. Subtracting the sky background 
fixes this. The sky background is calculated using the median pixel value in a 
coarse grid ``ntiles`` by ``ntiles`` on each image. The median-sampled coarse 
grid is then interpolated and resampled back to the resolution of the data array 
and subtracted from each image. This removes any strong gradients or vignetting 
that may dominate the shift calculation. 

Summary of object settings
~~~~~~~~~~~~~~~~~~~~~~~~~~

A summary of the object's settings can be seen using:

.. doctest-skip::

    >>> d.print_summary()


Computing per-image offsets
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With a ``Donuts`` object, the x/y translation required to align the
image into the reference image is given by the ``measure_shift`` method,
e.g.:

.. doctest-skip::

    >>> x, y = d.measure_shift(checkimage='image.fits')

x and y are returned in pixel units.
