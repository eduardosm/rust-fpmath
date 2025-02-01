use rand::Rng as _;

use super::{mkfloat, purify};
use crate::data::create_prng;

#[test]
fn test_round() {
    test_with(|arg| {
        let expected = fpmath::round(f64::from(arg));
        let actual = fpmath::round(arg);

        assert!(expected == f64::from(actual), "round({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_floor() {
    test_with(|arg| {
        let expected = fpmath::floor(f64::from(arg));
        let actual = fpmath::floor(arg);

        assert!(expected == f64::from(actual), "floor({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_ceil() {
    test_with(|arg| {
        let expected = fpmath::ceil(f64::from(arg));
        let actual = fpmath::ceil(arg);

        assert!(expected == f64::from(actual), "ceil({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_trunc() {
    test_with(|arg| {
        let expected = fpmath::trunc(f64::from(arg));
        let actual = fpmath::trunc(arg);

        assert!(expected == f64::from(actual), "trunc({arg:e}) = {actual:e}");
    });
}

fn test_with(test_f: fn(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        test_f(mkfloat(0, e, false));
        test_f(mkfloat(0, e, true));
        test_f(mkfloat(u32::MAX, e, false));
        test_f(mkfloat(u32::MAX, e, true));

        for _ in 0..10_000 {
            let m = rng.random::<u32>();
            test_f(mkfloat(m, e, true));
            test_f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=100_000 {
        let arg = arg as f32;
        test_f(purify(arg));
        test_f(purify(-arg));
        test_f(purify(arg + 0.25));
        test_f(purify(-arg + 0.25));
        test_f(purify(arg + 0.5));
        test_f(purify(-arg + 0.5));
        test_f(purify(arg + 0.75));
        test_f(purify(-arg + 0.75));
    }

    for e in 0..23 {
        for delta in -1000..=1000 {
            let arg = mkfloat(0, e, false) + delta as f32;
            test_f(purify(arg));
            test_f(purify(-arg));
            test_f(purify(arg + 0.25));
            test_f(purify(-arg + 0.25));
            test_f(purify(arg + 0.5));
            test_f(purify(-arg + 0.5));
            test_f(purify(arg + 0.75));
            test_f(purify(-arg + 0.75));
        }
    }
}
