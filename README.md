# fpmath

[![GitHub Actions Status](https://github.com/eduardosm/rust-fpmath/workflows/CI/badge.svg)](https://github.com/eduardosm/rust-fpmath/actions)
[![crates.io](https://img.shields.io/crates/v/fpmath.svg)](https://crates.io/crates/fpmath)
[![Documentation](https://docs.rs/fpmath/badge.svg)](https://docs.rs/fpmath)
![MSRV](https://img.shields.io/badge/rustc-1.70+-lightgray.svg)
[![License](https://img.shields.io/crates/l/fpmath.svg)](#license)

fpmath is a pure-Rust floating point library that implements math functions for
`f32` and `f64`.

## Features

* Less than 1 ULP error in most functions.
* `f32` versions do not internally use `f64` arithmetic.
* All functions are also implemented for soft-float single-precision and double precision.
* `no_std`
* All functions are tested for accuracy ([MPFR] is used to calculate reference values).
* The included `generator` crate can generate all magic constants used in the algorithms.

[MPFR]: https://www.mpfr.org/

### Included functions

* Absolute value
* Copy sign
* Rounding (to nearest, towards zero, towards infinity, towards negative
  infinity)
* Exponential in base e, 2 and 10
* Logarithm in base e, 2 and 10
* Trigonometric (sine, cosine, tangent) in radians, degrees and
  half-revolutions
* Inverse trigonometric (arcsine, arccosine, arctangent) in radians, degrees
  and half-revolutions
* Hyperbolic (sine, cosine, tangent)
* Inverse hyperbolic (arcsine, arccosine, arctangent)
* Square and cube root
* Hypotenuse
* Power (floating point and integer exponent)

## License

Licensed under either of

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  <http://www.apache.org/licenses/LICENSE-2.0>)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or
  <http://opensource.org/licenses/MIT>)

at your option.
