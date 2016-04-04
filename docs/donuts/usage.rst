******************
User Documentation
******************

To measure the shifts present in  a series of images requires the
definition of a reference frame. An image is chosen from the series
to be the reference. Then Donuts attempts to measure the offset of           
each subsequent image with respect to the reference frame. Reference
images should be free from unwanted artefacts (satellite trails etc.)

Initialization and Reference Image 
==================================

Initialization of the algorithm is done at the same time a reference 
image is generated. Below we import the Donuts class from the donuts
module and set up an instance with some suitable values.

```python
# import the Donuts class
from donuts import Donuts

# initialize the instance and generate a reference image
d=Donuts(refimage='reference.fits', image_ext=0, exposure='EXPTIME',
         normalise=True, subtract_bkg=True, prescan_width=0,
         overscan_width=0, border=64, ntiles=32)
``` 

To set up Donuts correctly you need to know the following information
about your imaging device. 

* ```image_ext``` : the fits image extension to use when measuring shifts. Default is 0. Ignore this option if your data contains only 1 fits extension.
* ```exposure``` : the header keyword decribing the exposure time of each image. Default is ```EXPTIME```. Ignore this option if ```normalise``` is set to ```FALSE```, see below.
* ```prescan_width``` : the width (in pixels) of your CCD's prescan region. Default is 0. Ignore this option if your CCD has no prescan region.
* ```overscan_width``` : the width (in pixels) of your CCD's overscan region. Default is 0. Ignore this option if your CCD has no overscan region. 
* ```border``` : the width of the exclusion border applied to the edges of the CCD. Border pixels typically have a non-linear response and are recommended to be excluded. Default is 64 pixels. 
* ```ntiles``` : the number of tiles used in the background subtraction, see below. Use more tiles if the sky background is varying rapidly. Default is 32. 

Donuts applies several corrections to each image before measuring the shifts. 

* 1. Donuts measures the shape of the image array. If an axis is found to be
indivisible by 16, that image axis is truncated to the nearest multiple of 512
pixels. This affects CCDs with non-standard image sizes only. 512 x 512 pixels 
is assumed to be the smallest usable CCD size.
* 2. Donuts then optionally (but highly recommended) removes a boarder from the 
edge of the image using the option described above.
* 3. An optional (but recommended) correction for exposure time is possible by 
setting the option ```normalise=True``` and providing the header keyword that 
contains the exposuretime, using the ```exposure``` option described above. 
This allows for varying exposure times without the need for separate analysis. 
* 4. A futher optional (but recommended) correction for the sky background is
possible. By design, Donuts will lock onto the most dominant 'feature' in a series
of images. Under normal circumstances this will be the forrest of stellar profiles
in the collapsed image projections. However, under new moon conditions for 
example, a strong vignetting pattern, might be the dominate feature and will 
result in poor performance. Subtracting the sky background fixes this. The sky
background is calculated using the median pixel value in a coarse grid
 ```ntiles``` by ```ntiles``` on each image. The median-sampled coarse grid is
then resampled back to the resolution of the data array and subtracted from each 
image. This removes any strong gradients or vignetting that may dominate the
shift calculation. 

Once the object has been initialized and a reference image has been generated, 
all subsequent images are corrected using the same process as the reference 
image. Finally the shift between each image and reference is measured. To measure 
the shift, simply pass the Donuts object a ```checkimage``` using the 
```measure_shift``` method. 

```python
d.measure_shift(checkimage='compare.fits')
```

The shift in X and Y is returned using dimmensions of ```astropy.units.pixel```.

A summary of the object's settings can be seen using:

```python
d.print_summary()
```
