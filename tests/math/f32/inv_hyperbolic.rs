use super::{check_result, mk_normal, purify};
use crate::create_prng;

#[test]
fn test_asinh() {
    let mut max_error: f32 = 0.0;
    test_asinh_with(|x| {
        let expected = fpmath::asinh(f64::from(x));
        let actual = fpmath::asinh(x);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_asinh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        f(mk_normal(0, e, false));
        f(mk_normal(0, e, true));
        f(mk_normal(super::MAX_MANTISSA, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, true));

        for _ in 0..20000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
            f(mk_normal(m, e, true));
        }
    }
}

#[test]
fn test_acosh() {
    let mut max_error: f32 = 0.0;
    test_acosh_with(|x| {
        let expected = fpmath::acosh(f64::from(x));
        let actual = fpmath::acosh(x);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_acosh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in 0..=127 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..40000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for e in -20..=-2 {
        for _ in 0..5000 {
            let m = super::gen_mantissa(&mut rng);
            let x = purify(1.0 + mk_normal(m, e, false));
            f(x);
        }
    }
}

#[test]
fn test_atanh() {
    let mut max_error: f32 = 0.0;
    test_atanh_with(|x| {
        let expected = fpmath::atanh(f64::from(x));
        let actual = fpmath::atanh(x);
        assert_total_eq!(fpmath::atanh(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atanh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    // Exhaustive test near ±1
    for m in 0..(1 << 23) {
        f(mk_normal(m, -1, false));
    }

    for e in -126..=-1 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..40000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }
}
