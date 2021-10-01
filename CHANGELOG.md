## [0.2.0] - 2021-10

### Added
- [`mappy`](https://pypi.org/project/mappy/) (Python bindings for minimap2) is used for parsing SAM/BAM files.
- 

### Fixed
- 

### Changed
- 

## [0.1.4] - 2021-01-05

### Added
- CI/CD support with GitHub Actions

### Fixed
- Segfault when loading alignments into Match objects

### Changed
- String attributes in Match instances are now char* types in C++ sam_module
- Travis CI is no longer the vendor for CI/CD

## [0.1.3] - 2021-01-03

### Added

- (#10) Merged the pull request that drastically improved the C++ extension by creating a Python-typed class object

### Changed
- Drastically reduced RAM by not instantiating objects for unaligned reads
