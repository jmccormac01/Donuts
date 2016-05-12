# Change log

## [Unreleased] - 12/05/2016

## Added

- Started [keeping a changelog]

## [0.0.1dev212] - 11/05/2016

### Changed

- Introduced the `Image` class, abstracting away the complex image data
  manipulation and cross correlation, leaving the `Donuts` interface the
  same. This has the advantage of allowing the user to access more "meta"
  information about the process.  Issue: [#23], PR: [#25]

## [0.0.1dev184] - 11/05/2016 

### Changed

- Split out the large test data files into an external repository, to
  ensure the `pypi` package is light. Issue: [#16], PR: [#18]

[Unreleased]: https://github.com/jmccormac01/Donuts/compare/282ca86d01ef...devel

[0.0.1dev212]: https://github.com/jmccormac01/Donuts/compare/4798806aa3ef...282ca86d01ef

[#23]: https://github.com/jmccormac01/Donuts/issues/23

[#25]: https://github.com/jmccormac01/Donuts/pull/25

[0.0.1dev184]: https://github.com/jmccormac01/Donuts/compare/94d68d5129c9...9d29cedebbfb

[#16]: https://github.com/jmccormac01/Donuts/issues/16

[#18]: https://github.com/jmccormac01/Donuts/pull/18

[keeping a changelog]: http://keepachangelog.com/
