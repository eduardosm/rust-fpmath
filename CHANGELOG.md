# Changelog

## 0.2.0 (unreleased)

### Breaking

- Soft-float support has been removed. Use the
  [`swfp` crate](https://crates.io/crates/swfp) instead.
- Some math functions have been renamed (`log` to `ln`, `log_1p` to `ln_1p`,
  `tgamma` to `gamma`, `lgamma` to `ln_gamma`)

### Fixed

- Fixed large `atanh` error near -1.
- Fixed large `ln_gamma` errors for some negative inputs.
- Symmetry of trigonometric functions now applies to the sign of zero too.
- Fixed `atanpi` and `atan2pi` exceeding 1 ULP of error.

### Changed

- Poles of `tand` and `tanpi` now have alternating signs.

## 0.1.1 (2024-10-14)

### Added

- Gamma (`tgamma`) and log-gamma (`lgamma`) functions.

## 0.1.0 (2024-09-17)

- Initial release
