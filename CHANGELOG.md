# Changelog

## 0.2.0 (unreleased)

### Breaking

- Soft-float support has been removed. Use the
  [`swfp` crate](https://crates.io/crates/swfp) instead.
- Some math functions have been renamed (`log` to `ln`, `log_1p` to `ln_1p`,
  `tgamma` to `gamma`, `lgamma` to `ln_gamma`)

### Fixed

- Fixed cases of errors being larger than documented maximum.
- Symmetry of trigonometric functions now applies to the sign of zero too.

### Changed

- Maximum error has been brought down to 0.55 ULP from 1 ULP.
- Poles of `tand` and `tanpi` now have alternating signs.
- `f32` math functions now use `f64` internally.

## 0.1.1 (2024-10-14)

### Added

- Gamma (`tgamma`) and log-gamma (`lgamma`) functions.

## 0.1.0 (2024-09-17)

- Initial release
