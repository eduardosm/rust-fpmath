use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, purify, select_threshold, RUG_PREC};
use crate::create_prng;

#[test]
fn test_asin_acos() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = rug::Float::with_val(RUG_PREC, x).asin();
        let expected_acos = rug::Float::with_val(RUG_PREC, x).acos();

        let actual_asin = fpmath::asin(x);
        let actual_acos = fpmath::acos(x);

        let asin_err = calc_error_ulp(actual_asin, expected_asin);
        let acos_err = calc_error_ulp(actual_acos, expected_acos);

        max_asin_error = max_asin_error.max(asin_err);
        max_acos_error = max_acos_error.max(acos_err);

        let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
        assert!(
            asin_err < asin_threshold,
            "asin({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
        );

        let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
        assert!(
            acos_err < acos_threshold,
            "acos({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
        );
    });
    eprintln!("max asin error = {max_asin_error}");
    eprintln!("max acos error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

#[test]
fn test_asind_acosd() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = rug::Float::with_val(RUG_PREC, x).asin_u(360);
        let expected_acos = rug::Float::with_val(RUG_PREC, x).acos_u(360);

        let actual_asin = fpmath::asind(x);
        let actual_acos = fpmath::acosd(x);

        let asin_err = calc_error_ulp(actual_asin, expected_asin);
        let acos_err = calc_error_ulp(actual_acos, expected_acos);

        max_asin_error = max_asin_error.max(asin_err);
        max_acos_error = max_acos_error.max(acos_err);

        let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
        assert!(
            asin_err < asin_threshold,
            "asind({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
        );

        let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
        assert!(
            acos_err < acos_threshold,
            "acosd({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
        );
    });
    eprintln!("max asind error = {max_asin_error}");
    eprintln!("max acosd error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

#[test]
fn test_asinpi_acospi() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = rug::Float::with_val(RUG_PREC, x).asin_pi();
        let expected_acos = rug::Float::with_val(RUG_PREC, x).acos_pi();

        let actual_asin = fpmath::asinpi(x);
        let actual_acos = fpmath::acospi(x);

        let asin_err = calc_error_ulp(actual_asin, expected_asin);
        let acos_err = calc_error_ulp(actual_acos, expected_acos);

        max_asin_error = max_asin_error.max(asin_err);
        max_acos_error = max_acos_error.max(acos_err);

        let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
        assert!(
            asin_err < asin_threshold,
            "asinpi({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
        );

        let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
        assert!(
            acos_err < acos_threshold,
            "acospi({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
        );
    });
    eprintln!("max asinpi error = {max_asin_error}");
    eprintln!("max acospi error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

fn test_asin_acos_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=-1 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u64::MAX, e, false));
        f(mkfloat(u64::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u64>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }

    for i in 1..=1000 {
        let x = (i as f64) / 1000.0;
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
        assert_eq!(purify(fpmath::atan(-x)), purify(-actual));

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atan error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atand() {
    let mut max_error: f64 = 0.0;
    test_atan_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_u(360);
        let actual = fpmath::atand(x);
        assert_eq!(purify(fpmath::atand(-x)), purify(-actual));

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atand({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atand error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atanpi() {
    let mut max_error: f64 = 0.0;
    test_atan_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_pi();
        let actual = fpmath::atanpi(x);
        assert_eq!(purify(fpmath::atanpi(-x)), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atanpi({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atanpi error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atan_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
}

#[test]
fn test_atan2() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected = rug::Float::with_val(RUG_PREC, y).atan2(&rug::Float::with_val(RUG_PREC, x));
        let actual = fpmath::atan2(y, x);
        assert_eq!(purify(fpmath::atan2(-y, x)), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atan2({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2d() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_u(&rug::Float::with_val(RUG_PREC, x), 360);
        let actual = fpmath::atan2d(y, x);
        assert_eq!(purify(fpmath::atan2d(-y, x)), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atan2d({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2d error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2pi() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_pi(&rug::Float::with_val(RUG_PREC, x));
        let actual = fpmath::atan2pi(y, x);
        assert_eq!(purify(fpmath::atan2pi(-y, x)), purify(-actual));

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.95, 1.9);
        assert!(
            err < threshold,
            "atan2pi({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2pi error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atan2_with(mut f: impl FnMut(f64, f64)) {
    let mut rng = create_prng();

    for ey in -1022..=1023 {
        if matches!(ey, -900..=900) && (ey & 3) != 3 {
            continue; // speed up tests
        }
        for ex in -1022..=1023 {
            if matches!(ex, -900..=900) && (ex & 3) != 3 {
                continue; // speed up tests
            }
            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            f(mkfloat(my, ey, false), mkfloat(mx, ex, false));
            f(mkfloat(my, ey, false), mkfloat(mx, ex, true));
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, e, false), mkfloat(mx, e, sx));

            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, 0, false), mkfloat(mx, e, sx));

            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, e, false), mkfloat(mx, 0, sx));
        }
    }
}
