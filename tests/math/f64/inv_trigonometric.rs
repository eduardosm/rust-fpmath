use super::{RUG_PREC, calc_error_ulp, mk_normal, purify};
use crate::create_prng;

#[test]
fn test_asin() {
    let mut max_error: f64 = 0.0;
    test_asin_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin();
        let actual = fpmath::asin(x);
        assert_total_eq!(fpmath::asin(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "asin({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_acos() {
    let mut max_error: f64 = 0.0;
    test_acos_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos();
        let actual = fpmath::acos(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "acos({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_asind() {
    let mut max_error: f64 = 0.0;
    test_asin_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin_u(360);
        let actual = fpmath::asind(x);
        assert_total_eq!(fpmath::asind(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "asind({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_acosd() {
    let mut max_error: f64 = 0.0;
    test_acos_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos_u(360);
        let actual = fpmath::acosd(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "acosd({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.4999);
}

#[test]
fn test_asinpi() {
    let mut max_error: f64 = 0.0;
    test_asin_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin_pi();
        let actual = fpmath::asinpi(x);
        assert_total_eq!(fpmath::asinpi(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "asinpi({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_acospi() {
    let mut max_error: f64 = 0.0;
    test_acos_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos_pi();
        let actual = fpmath::acospi(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "acospi({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error >= 0.5);
}

fn test_asin_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=-1 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..10000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for i in 1..=1000 {
        let x = purify((i as f64) / 1000.0);
        f(x);
    }
}

fn test_acos_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=-1 {
        f(mk_normal(0, e, false));
        f(mk_normal(0, e, true));
        f(mk_normal(super::MAX_MANTISSA, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, true));

        for _ in 0..10000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
            f(mk_normal(m, e, true));
        }
    }

    for i in 1..=1000 {
        let x = purify((i as f64) / 1000.0);
        f(x);
        f(-x);
    }
}

#[test]
fn test_atan() {
    let mut max_error: f64 = 0.0;
    test_atan_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan();
        let actual = fpmath::atan(x);
        assert_total_eq!(fpmath::atan(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "atan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atand() {
    let mut max_error: f64 = 0.0;
    test_atan_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_u(360);
        let actual = fpmath::atand(x);
        assert_total_eq!(fpmath::atand(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "atand({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atanpi() {
    let mut max_error: f64 = 0.0;
    test_atan_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_pi();
        let actual = fpmath::atanpi(x);
        assert_total_eq!(fpmath::atanpi(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "atanpi({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atan_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..20000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }
}

#[test]
fn test_atan2() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected = rug::Float::with_val(RUG_PREC, y).atan2(&rug::Float::with_val(RUG_PREC, x));
        let actual = fpmath::atan2(y, x);
        assert_total_eq!(fpmath::atan2(-y, x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "atan2({y:e}, {x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2d() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_u(&rug::Float::with_val(RUG_PREC, x), 360);
        let actual = fpmath::atan2d(y, x);
        assert_total_eq!(fpmath::atan2d(-y, x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "atan2d({y:e}, {x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2pi() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_pi(&rug::Float::with_val(RUG_PREC, x));
        let actual = fpmath::atan2pi(y, x);
        assert_total_eq!(fpmath::atan2pi(-y, x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "atan2pi({y:e}, {x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atan2_with(mut f: impl FnMut(f64, f64)) {
    let mut rng = create_prng();

    f(f64::INFINITY, f64::INFINITY);
    f(f64::INFINITY, f64::NEG_INFINITY);

    for ey in -1022..=1023 {
        if matches!(ey, -900..=900) && (ey & 3) != 3 {
            continue; // speed up tests
        }
        for ex in -1022..=1023 {
            if matches!(ex, -900..=900) && (ex & 3) != 3 {
                continue; // speed up tests
            }
            let my = super::gen_mantissa(&mut rng);
            let mx = super::gen_mantissa(&mut rng);
            f(mk_normal(my, ey, false), mk_normal(mx, ex, false));
            f(mk_normal(my, ey, false), mk_normal(mx, ex, true));
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let my = super::gen_mantissa(&mut rng);
            let mx = super::gen_mantissa(&mut rng);
            f(mk_normal(my, e, false), mk_normal(mx, e, false));
            f(mk_normal(my, e, false), mk_normal(mx, e, true));

            let my = super::gen_mantissa(&mut rng);
            let mx = super::gen_mantissa(&mut rng);
            f(mk_normal(my, 0, false), mk_normal(mx, e, false));
            f(mk_normal(my, 0, false), mk_normal(mx, e, true));

            let my = super::gen_mantissa(&mut rng);
            let mx = super::gen_mantissa(&mut rng);
            f(mk_normal(my, e, false), mk_normal(mx, 0, false));
            f(mk_normal(my, e, false), mk_normal(mx, 0, true));
        }
    }
}
