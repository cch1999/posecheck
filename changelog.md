# Changelog

## [1.1] - 18-01-2024

### Changed
- [IMPORTANT - PLEASE READ IF YOU INTEND TO PUBLISH] Changed the way strain energy is calculated. Previously, the strain energy was calculated as the difference between the energy of the raw output and a fully relaxed structure. This meant that the energy term of the raw output exploded due to slight imperfections in local bond lengths and angles and did not reflect global strain well. Now we perform a small bit of local relaxation (max displacement in the atom positions of 0.1 Å) and use the energy of this structure as the reference energy. This means that the strain energy is now a more accurate reflection of the global strain of the molecule.

### Added
- Data and notebooks needed to reproduce plots from the paper.
- Validator that warms the user if `reduce` cannot be found.
- CI tests
- Full environment specs

### Fixed
- Improved doc string coverage

## [1.0] - 24-09-2023

Initial release of the software with basic features for chemical analysis.
