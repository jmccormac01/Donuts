## Donuts
[https://travis-ci.org/jmccormac01/Donuts.svg?branch=master]()

A science frame autoguiding and image alignment algorithm with sub-pixel precision, capable of guiding on defocused stars. 

## Motivation

We operate and have access to several telescopes (NGTS, NITES, Warwick 1m, 1.5m San Pedro Martir) that require precise autoguiding. Sometimes we must defocus the telescope for various reasons, but we'd still like to autoguide. Donuts was designed to allow this. It had to be simple, fast and accurate. It works well as an autoguiding algorithm for equatorial telescopes (no field rotation). 

The process for aligning apertures for photometry is essentially the reverse. Rather than correcting the telescope pointing, the apertures must track the drift of the stars. Donuts can therefore be used to track the stellar positions for CCD photometry also. 

The algorithm has its limitations. It currently does not deal with rotation and large drifts, where the field moves by approx. half a FOV or more. Our paper describing the details can be found here:

http://adsabs.harvard.edu/abs/2013PASP..125..548M

## Code Example

 
## Installation

This package will hopefully be an Astropy affiliated package when I am done

## API Reference

Documentation will be done with Sphinx - not complete yet.

## Tests

There is currently a very simple test built into main. This takes three images and measures the shift between them, using the first as the reference. The expected results are:

|        *Image ID*           | X shift (pix) |  Y shift (pix) |
| --------------------------- | ------------- | -------------- |
| IMAGE80520160114005520.fits |    -0.09      |      0.24      |
| IMAGE80520160114005533.fits |     0.01      |      0.14      |

The corrections, are those required to adjust the stars back to their location in the reference image.

## Contributors

[@jmccormac01](https://github.com/jmccormac01) & [@mindriot101](https://github.com/mindriot101)


## License

MIT License

Copyright (c) 2016 James McCormac & Simon Walker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


