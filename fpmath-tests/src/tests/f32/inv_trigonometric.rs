use rand::Rng as _;

use super::{mkfloat, purify, select_threshold, RefResult};
use crate::data::create_prng;

#[test]
fn test_asin_acos() {
    let mut max_asin_error: f32 = 0.0;
    let mut max_acos_error: f32 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = RefResult::from_f64(fpmath::asin(f64::from(x)));
        let expected_acos = RefResult::from_f64(fpmath::acos(f64::from(x)));

        let actual_asin = fpmath::asin(x);
        let actual_acos = fpmath::acos(x);

        let asin_err = expected_asin.calc_error(actual_asin);
        let acos_err = expected_acos.calc_error(actual_acos);

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
    let mut max_asin_error: f32 = 0.0;
    let mut max_acos_error: f32 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = RefResult::from_f64(fpmath::asind(f64::from(x)));
        let expected_acos = RefResult::from_f64(fpmath::acosd(f64::from(x)));

        let actual_asin = fpmath::asind(x);
        let actual_acos = fpmath::acosd(x);

        let asin_err = expected_asin.calc_error(actual_asin);
        let acos_err = expected_acos.calc_error(actual_acos);

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
    let mut max_asin_error: f32 = 0.0;
    let mut max_acos_error: f32 = 0.0;
    test_asin_acos_with(|x| {
        let expected_asin = RefResult::from_f64(fpmath::asinpi(f64::from(x)));
        let expected_acos = RefResult::from_f64(fpmath::acospi(f64::from(x)));

        let actual_asin = fpmath::asinpi(x);
        let actual_acos = fpmath::acospi(x);

        let asin_err = expected_asin.calc_error(actual_asin);
        let acos_err = expected_acos.calc_error(actual_acos);

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

fn test_asin_acos_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=-1 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u32::MAX, e, false));
        f(mkfloat(u32::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }

    for i in 1..=1000 {
        let x = (i as f32) / 1000.0;
        f(x);
        f(-x);
    }
}

#[test]
fn test_atan() {
    let mut max_error: f32 = 0.0;
    test_atan_with(|x| {
        let expected = RefResult::from_f64(fpmath::atan(f64::from(x)));
        let actual = fpmath::atan(x);
        assert_eq!(purify(fpmath::atan(-x)), purify(-actual));

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atan error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atand() {
    let mut max_error: f32 = 0.0;
    test_atan_with(|x| {
        let expected = RefResult::from_f64(fpmath::atand(f64::from(x)));
        let actual = fpmath::atand(x);
        assert_eq!(purify(fpmath::atand(-x)), purify(-actual));

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atand({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atand error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atanpi() {
    let mut max_error: f32 = 0.0;
    test_atan_with(|x| {
        let expected = RefResult::from_f64(fpmath::atanpi(f64::from(x)));
        let actual = fpmath::atanpi(x);
        assert_eq!(purify(fpmath::atanpi(-x)), purify(-actual));

        let err = expected.calc_error(actual);
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

fn test_atan_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(u32::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
        }
    }
}

#[test]
fn test_atan2() {
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = RefResult::from_f64(fpmath::atan2(f64::from(y), f64::from(x)));
        let actual = fpmath::atan2(y, x);
        assert_eq!(purify(fpmath::atan2(-y, x)), purify(-actual));

        let err = expected.calc_error(actual);
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
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = RefResult::from_f64(fpmath::atan2d(f64::from(y), f64::from(x)));
        let actual = fpmath::atan2d(y, x);
        assert_eq!(purify(fpmath::atan2d(-y, x)), purify(-actual));

        let err = expected.calc_error(actual);
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
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = RefResult::from_f64(fpmath::atan2pi(f64::from(y), f64::from(x)));
        let actual = fpmath::atan2pi(y, x);
        assert_eq!(purify(fpmath::atan2pi(-y, x)), purify(-actual));

        let err = expected.calc_error(actual);
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

fn test_atan2_with(mut f: impl FnMut(f32, f32)) {
    let mut rng = create_prng();

    for ey in -126..=127 {
        for ex in -126..=127 {
            let my = rng.random::<u32>();
            let mx = rng.random::<u32>();
            f(mkfloat(my, ey, false), mkfloat(mx, ex, false));
            f(mkfloat(my, ey, false), mkfloat(mx, ex, true));
        }
    }

    for e in -126..=127 {
        for _ in 0..5000 {
            let my = rng.random::<u32>();
            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, e, false), mkfloat(mx, e, sx));

            let my = rng.random::<u32>();
            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, 0, false), mkfloat(mx, e, sx));

            let my = rng.random::<u32>();
            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            f(mkfloat(my, e, false), mkfloat(mx, 0, sx));
        }
    }
}
