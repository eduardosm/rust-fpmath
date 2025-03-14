use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, RUG_PREC};
use crate::create_prng;

#[test]
fn test_log() {
    let mut max_error: f64 = 0.0;
    test_log_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).ln();
        let actual = fpmath::log(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log_1p() {
    let mut max_error: f64 = 0.0;
    test_log1p_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).ln_1p();
        let actual = fpmath::log_1p(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log_1p({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log_1p error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log2() {
    let mut max_error: f64 = 0.0;
    test_log_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).log2();
        let actual = fpmath::log2(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log2({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log10() {
    let mut max_error: f64 = 0.0;
    test_log_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).log10();
        let actual = fpmath::log10(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log10({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log10 error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_log_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -100..=100 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..9000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..1000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f64);
    }

    f(f64::MIN_POSITIVE);
    f(f64::MAX);

    // subnormals
    for i in 0..52 {
        f(f64::from_bits(1 << i));
        f(f64::from_bits((1 << (i + 1)) - 1));
    }
}

fn test_log1p_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -100..=100 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..9000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..1000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f64);
    }

    f(f64::MIN_POSITIVE);
    f(f64::MAX);

    // 1 < x < 0
    for e in -1022..=-1 {
        f(mkfloat(0, e, true));
        f(mkfloat(u64::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, true));
        }
    }

    // subnormals
    for i in 0..52 {
        f(f64::from_bits(1 << i));
        f(-f64::from_bits(1 << i));
        f(f64::from_bits((1 << (i + 1)) - 1));
        f(-f64::from_bits((1 << (i + 1)) - 1));
    }
}
