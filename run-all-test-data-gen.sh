#!/bin/sh
set -e

DATA_DIR="$(dirname "$(readlink -f "$0")")/fpmath-tests/data"
if [ -e "$DATA_DIR" ]; then
    rm -r "$DATA_DIR"
fi

./run-test-data-gen.sh \
    f32::acosh \
    f32::asin_acos \
    f32::asind_acosd \
    f32::asinh \
    f32::asinpi_acospi \
    f32::atan \
    f32::atan2 \
    f32::atan2d \
    f32::atan2pi \
    f32::atand \
    f32::atanh \
    f32::atanpi \
    f32::cbrt \
    f32::exp \
    f32::exp10 \
    f32::exp2 \
    f32::hypot \
    f32::log \
    f32::log10 \
    f32::log2 \
    f32::log_1p \
    f32::pow \
    f32::powi \
    f32::sin_cos \
    f32::sind_cosd \
    f32::sinh_cosh \
    f32::sinpi_cospi \
    f32::sqrt \
    f32::tan \
    f32::tand \
    f32::tanh \
    f32::tanpi \
    f64::acosh \
    f64::asin_acos \
    f64::asind_acosd \
    f64::asinh \
    f64::asinpi_acospi \
    f64::atan \
    f64::atan2 \
    f64::atan2d \
    f64::atan2pi \
    f64::atand \
    f64::atanh \
    f64::atanpi \
    f64::cbrt \
    f64::exp \
    f64::exp10 \
    f64::exp2 \
    f64::hypot \
    f64::log \
    f64::log10 \
    f64::log2 \
    f64::log_1p \
    f64::pow \
    f64::powi \
    f64::sin_cos \
    f64::sind_cosd \
    f64::sinh_cosh \
    f64::sinpi_cospi \
    f64::sqrt \
    f64::tan \
    f64::tand \
    f64::tanh \
    f64::tanpi
