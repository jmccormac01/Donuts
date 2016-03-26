## Donuts

A science frame autoguiding and image alignment algorithm with sub-pixel precision, capable of guiding on defocused stars. 

## Motivation

We operate and have access to several telescopes (NGTS, NITES, Warwick 1m, 1.5m San Pedro Martir) that require precise autoguiding. Sometimes we must defocus the telescope for various reasons, but we'd still like to autoguide. Donuts was designed to allow this. It had to be simple, fast and accurate. It works well as an autoguiding algorithm for equatorial telescopes (no field rotation). 

The process for aligning apertures for photometry is essentially the reverse. Rather than the telescope pointing being corrected. The apertures must track the drift of the stars. Donuts can therefore be used to track the stellar positions for CCD photometry also. 

The algorithm has its limitations. It currently does not deal with rotation and large drifts, where the field moves by approx. half a FOV or more. Our paper describing the details can be found here:

http://adsabs.harvard.edu/abs/2013PASP..125..548M

## Code Example

 
## Installation

This package will hopefully be an Astropy affiliated package when I am done

## API Reference

Documentation will be done with Sphinx - not complete yet.

## Tests

There is currently a very simple test built into main. This takes three images and measures the shift between them, using the first as the reference. The expected results are:

| IMAGE80520160114005520.fits | X: -0.09 | Y: 0.24 |
| IMAGE80520160114005533.fits | X: 0.01  | Y: 0.14 |

The corrections, are those required to adjust the stars back to their location in the reference image.

## Contributors

@jmccormac01 & @mindriot101

## License

BSD
)
