use super::{RUG_PREC, check_result, mk_normal, mk_subnormal, purify};
use crate::create_prng;

#[test]
fn test_sin_cos() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) =
            rug::Float::with_val(RUG_PREC, x).sin_cos(rug::Float::new(RUG_PREC));

        let actual_sin1 = fpmath::sin(x);
        let actual_cos1 = fpmath::cos(x);
        let (actual_sin2, actual_cos2) = fpmath::sin_cos(x);
        assert_total_eq!(fpmath::sin(-x), -actual_sin1);
        assert_total_eq!(fpmath::cos(-x), actual_cos1);
        assert_total_eq!(fpmath::sin_cos(-x), (-actual_sin2, actual_cos2));

        check_result(
            x,
            actual_sin1,
            expected_sin.clone(),
            0.51,
            &mut max_sin1_error,
        );
        check_result(x, actual_sin2, expected_sin, 0.51, &mut max_sin2_error);
        check_result(
            x,
            actual_cos1,
            expected_cos.clone(),
            0.51,
            &mut max_cos1_error,
        );
        check_result(x, actual_cos2, expected_cos, 0.51, &mut max_cos2_error);
    });
    eprintln!("max sin1 error = {max_sin1_error}");
    eprintln!("max sin2 error = {max_sin2_error}");
    eprintln!("max cos1 error = {max_cos1_error}");
    eprintln!("max cos2 error = {max_cos2_error}");
    assert!(max_sin1_error >= 0.5);
    assert!(max_sin2_error >= 0.5);
    assert!(max_cos1_error >= 0.5);
    assert!(max_cos2_error >= 0.5);
}

#[test]
fn test_sind_cosd() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    test_with(|x| {
        let expected_sin = rug::Float::with_val(RUG_PREC, x).sin_u(360);
        let expected_cos = rug::Float::with_val(RUG_PREC, x).cos_u(360);

        let actual_sin1 = fpmath::sind(x);
        let actual_cos1 = fpmath::cosd(x);
        let (actual_sin2, actual_cos2) = fpmath::sind_cosd(x);
        assert_total_eq!(fpmath::sind(-x), -actual_sin1);
        assert_total_eq!(fpmath::cosd(-x), actual_cos1);
        assert_total_eq!(fpmath::sind_cosd(-x), (-actual_sin2, actual_cos2));

        check_result(
            x,
            actual_sin1,
            expected_sin.clone(),
            0.51,
            &mut max_sin1_error,
        );
        check_result(x, actual_sin2, expected_sin, 0.51, &mut max_sin2_error);
        check_result(
            x,
            actual_cos1,
            expected_cos.clone(),
            0.51,
            &mut max_cos1_error,
        );
        check_result(x, actual_cos2, expected_cos, 0.51, &mut max_cos2_error);
    });
    eprintln!("max sind1 error = {max_sin1_error}");
    eprintln!("max sind2 error = {max_sin2_error}");
    eprintln!("max cosd1 error = {max_cos1_error}");
    eprintln!("max cosd2 error = {max_cos2_error}");
    assert!(max_sin1_error >= 0.5);
    assert!(max_sin2_error >= 0.5);
    assert!(max_cos1_error >= 0.5);
    assert!(max_cos2_error >= 0.5);
}

#[test]
fn test_sinpi_cospi() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    test_with(|x| {
        let expected_sin = rug::Float::with_val(RUG_PREC, x).sin_pi();
        let expected_cos = rug::Float::with_val(RUG_PREC, x).cos_pi();

        let actual_sin1 = fpmath::sinpi(x);
        let actual_cos1 = fpmath::cospi(x);
        let (actual_sin2, actual_cos2) = fpmath::sinpi_cospi(x);
        assert_total_eq!(fpmath::sinpi(-x), -actual_sin1);
        assert_total_eq!(fpmath::cospi(-x), actual_cos1);
        assert_total_eq!(fpmath::sinpi_cospi(-x), (-actual_sin2, actual_cos2));

        check_result(
            x,
            actual_sin1,
            expected_sin.clone(),
            0.51,
            &mut max_sin1_error,
        );
        check_result(x, actual_sin2, expected_sin, 0.51, &mut max_sin2_error);
        check_result(
            x,
            actual_cos1,
            expected_cos.clone(),
            0.51,
            &mut max_cos1_error,
        );
        check_result(x, actual_cos2, expected_cos, 0.51, &mut max_cos2_error);
    });
    eprintln!("max sinpi1 error = {max_sin1_error}");
    eprintln!("max sinpi2 error = {max_sin2_error}");
    eprintln!("max cospi1 error = {max_cos1_error}");
    eprintln!("max cospi2 error = {max_cos2_error}");
    assert!(max_sin1_error >= 0.5);
    assert!(max_sin2_error >= 0.5);
    assert!(max_cos1_error >= 0.5);
    assert!(max_cos2_error >= 0.5);
}

#[test]
fn test_tan() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).tan();
        let actual = fpmath::tan(x);
        assert_total_eq!(fpmath::tan(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tand() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).tan_u(360);
        let actual = fpmath::tand(x);
        assert_total_eq!(fpmath::tand(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tanpi() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).tan_pi();
        let actual = fpmath::tanpi(x);
        assert_total_eq!(fpmath::tanpi(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for _ in 0..100_000 {
        let m = super::gen_mantissa(&mut rng);
        f(mk_subnormal(m, false));
    }

    for e in -1022..=1023 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..5000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for arg in 1..=200_000 {
        f(arg as f64);
    }

    // Problematic value in
    // "ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit"
    f(1.0e22);

    // Magic values from rust-libm tests
    f(3.141592025756836);
    f(3.141592033207416);
    f(3.141592144966125);
    f(3.141592979431152);
    f(3054214.5490637687);
    f(917340800458.2274);

    let f2s = [1.0, (1 << 20) as f64];
    let f3s = [
        1.0, 1.1, 1.01, 1.001, 1.0001, 1.00001, 1.000001, 1.0000001, 0.9, 0.99, 0.999, 0.9999,
        0.99999, 0.999999, 0.9999999,
    ];

    for f1 in 1..=100 {
        for f2 in f2s {
            for f3 in f3s {
                f(purify(std::f64::consts::FRAC_PI_8 * (f1 as f64) * f2 * f3));
            }
        }
    }
}
