use super::{calc_error_ulp, mk_normal, mk_subnormal, purify};
use crate::create_prng;

#[test]
fn test_sin_cos() {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) = fpmath::sin_cos(f64::from(x));

        let actual_sin1 = fpmath::sin(x);
        let actual_cos1 = fpmath::cos(x);
        let (actual_sin2, actual_cos2) = fpmath::sin_cos(x);
        assert_total_eq!(fpmath::sin(-x), -actual_sin1);
        assert_total_eq!(fpmath::cos(-x), actual_cos1);
        assert_total_eq!(fpmath::sin_cos(-x), (-actual_sin2, actual_cos2));

        let sin1_err = calc_error_ulp(actual_sin1, expected_sin);
        let sin2_err = calc_error_ulp(actual_sin2, expected_sin);
        let cos1_err = calc_error_ulp(actual_cos1, expected_cos);
        let cos2_err = calc_error_ulp(actual_cos2, expected_cos);

        max_sin1_error = max_sin1_error.max(sin1_err);
        max_sin2_error = max_sin2_error.max(sin2_err);
        max_cos1_error = max_cos1_error.max(cos1_err);
        max_cos2_error = max_cos2_error.max(cos2_err);

        assert!(
            sin1_err < 0.51,
            "sin({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
        );

        assert!(
            sin2_err < 0.51,
            "sin({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
        );

        assert!(
            cos1_err < 0.51,
            "cos({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
        );

        assert!(
            cos2_err < 0.51,
            "cos({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
        );
    });
    eprintln!("max sin1 error = {max_sin1_error}");
    eprintln!("max sin2 error = {max_sin2_error}");
    eprintln!("max cos1 error = {max_cos1_error}");
    eprintln!("max cos2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_sind_cosd() {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) = fpmath::sind_cosd(f64::from(x));

        let actual_sin1 = fpmath::sind(x);
        let actual_cos1 = fpmath::cosd(x);
        let (actual_sin2, actual_cos2) = fpmath::sind_cosd(x);
        assert_total_eq!(fpmath::sind(-x), -actual_sin1);
        assert_total_eq!(fpmath::cosd(-x), actual_cos1);
        assert_total_eq!(fpmath::sind_cosd(-x), (-actual_sin2, actual_cos2));

        let sin1_err = calc_error_ulp(actual_sin1, expected_sin);
        let sin2_err = calc_error_ulp(actual_sin2, expected_sin);
        let cos1_err = calc_error_ulp(actual_cos1, expected_cos);
        let cos2_err = calc_error_ulp(actual_cos2, expected_cos);

        max_sin1_error = max_sin1_error.max(sin1_err);
        max_sin2_error = max_sin2_error.max(sin2_err);
        max_cos1_error = max_cos1_error.max(cos1_err);
        max_cos2_error = max_cos2_error.max(cos2_err);

        assert!(
            sin1_err < 0.51,
            "sind({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
        );
        assert!(
            sin2_err < 0.51,
            "sind({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
        );

        assert!(
            cos1_err < 0.51,
            "cosd({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
        );
        assert!(
            cos2_err < 0.51,
            "cosd({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
        );
    });
    eprintln!("max sind1 error = {max_sin1_error}");
    eprintln!("max sind2 error = {max_sin2_error}");
    eprintln!("max cosd1 error = {max_cos1_error}");
    eprintln!("max cosd2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_sinpi_cospi() {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) = fpmath::sinpi_cospi(f64::from(x));

        let actual_sin1 = fpmath::sinpi(x);
        let actual_cos1 = fpmath::cospi(x);
        let (actual_sin2, actual_cos2) = fpmath::sinpi_cospi(x);
        assert_total_eq!(fpmath::sinpi(-x), -actual_sin1);
        assert_total_eq!(fpmath::cospi(-x), actual_cos1);
        assert_total_eq!(fpmath::sinpi_cospi(-x), (-actual_sin2, actual_cos2));

        let sin1_err = calc_error_ulp(actual_sin1, expected_sin);
        let sin2_err = calc_error_ulp(actual_sin2, expected_sin);
        let cos1_err = calc_error_ulp(actual_cos1, expected_cos);
        let cos2_err = calc_error_ulp(actual_cos2, expected_cos);

        max_sin1_error = max_sin1_error.max(sin1_err);
        max_sin2_error = max_sin2_error.max(sin2_err);
        max_cos1_error = max_cos1_error.max(cos1_err);
        max_cos2_error = max_cos2_error.max(cos2_err);

        assert!(
            sin1_err < 0.51,
            "sinpi({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
        );
        assert!(
            sin2_err < 0.51,
            "sinpi({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
        );

        assert!(
            cos1_err < 0.51,
            "cospi({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
        );
        assert!(
            cos2_err < 0.51,
            "cospi({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
        );
    });
    eprintln!("max sinpi1 error = {max_sin1_error}");
    eprintln!("max sinpi2 error = {max_sin2_error}");
    eprintln!("max cospi1 error = {max_cos1_error}");
    eprintln!("max cospi2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_tan() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tan(f64::from(x));
        let actual = fpmath::tan(x);
        assert_total_eq!(fpmath::tan(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "tan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tand() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tand(f64::from(x));
        let actual = fpmath::tand(x);
        assert_total_eq!(fpmath::tand(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "tand({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tanpi() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tanpi(f64::from(x));
        let actual = fpmath::tanpi(x);
        assert_total_eq!(fpmath::tanpi(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "tanpi({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for m in 0..(1 << 23) {
        f(mk_subnormal(m, false));
    }

    for e in -126..=127 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..50000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for arg in 1..=20_000 {
        f(arg as f32);
    }

    // Problematic value in
    // "ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit"
    f(1.0e22);

    let f2s = [1.0, (1 << 9) as f32];
    let f3s = [
        1.0, 1.1, 1.01, 1.001, 1.0001, 1.00001, 1.000001, 1.0000001, 0.9, 0.99, 0.999, 0.9999,
        0.99999, 0.999999, 0.9999999,
    ];

    for f1 in 1..=100 {
        for f2 in f2s {
            for f3 in f3s {
                f(purify(std::f32::consts::FRAC_PI_8 * (f1 as f32) * f2 * f3));
            }
        }
    }
}
