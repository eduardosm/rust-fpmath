use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, purify, purify2};
use crate::data::create_prng;

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
        assert_eq!(purify(fpmath::sin(-x)), purify(-actual_sin1));
        assert_eq!(purify(fpmath::cos(-x)), purify(actual_cos1));
        assert_eq!(
            purify2(fpmath::sin_cos(-x)),
            purify2((-actual_sin2, actual_cos2))
        );

        let sin1_err = calc_error_ulp(actual_sin1, expected_sin);
        let sin2_err = calc_error_ulp(actual_sin2, expected_sin);
        let cos1_err = calc_error_ulp(actual_cos1, expected_cos);
        let cos2_err = calc_error_ulp(actual_cos2, expected_cos);

        max_sin1_error = max_sin1_error.max(sin1_err);
        max_sin2_error = max_sin2_error.max(sin2_err);
        max_cos1_error = max_cos1_error.max(cos1_err);
        max_cos2_error = max_cos2_error.max(cos2_err);

        assert!(
            sin1_err < 0.9,
            "sin({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
        );

        assert!(
            sin2_err < 0.9,
            "sin({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
        );

        assert!(
            cos1_err < 0.9,
            "cos({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
        );

        assert!(
            cos2_err < 0.9,
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
fn test_tan() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tan(f64::from(x));
        let actual = fpmath::tan(x);
        assert_eq!(purify(fpmath::tan(-x)), purify(-actual));

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "tan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max tan error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(u32::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=1000 {
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
