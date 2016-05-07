********
Advanced
********

For those who want more control over how the data is transformed, we provide a
number of places to hook into the process.

``Image`` class
---------------

This class really encapsulates the process of transforming a science image, and
computing the separation between processed images. The ``Donuts`` class is merely
a high level API over the ``Image`` class.

The easiest way to customise the transformation process (trimming, normalisation
and background subtraction) is to subclass ``Image`` and override the
``preconstruct_hook`` and ``postconstruct_hook`` methods, for example::

    from donuts.image import Image

    class MyCustomImage(Image):
        def preconstruct_hook(self):
            self.med_image = np.median(self.raw_image)
            self.raw_image -= self.med_image

        def postconstruct_hook(self):
            self.backsub += self.med_image

Then to "install" this custom class into the process, construct a ``Donuts``
object using the ``image_class`` parameter of the constructor::

    from donuts import Donuts

    d = Donuts(filename, image_class=MyCustomImage)


If further customisation is required, override other methods in the process but
currently no convenience hooks are present to make the process easier. As always
feel free to check out the source code to see what's going on.
