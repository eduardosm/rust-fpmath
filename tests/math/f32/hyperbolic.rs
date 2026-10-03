use super::{check_result, mk_normal};
use crate::create_prng;

#[test]
fn test_sinh_cosh() {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) = fpmath::sinh_cosh(f64::from(x));

        let actual_sin1 = fpmath::sinh(x);
        let actual_cos1 = fpmath::cosh(x);
        let (actual_sin2, actual_cos2) = fpmath::sinh_cosh(x);
        assert_total_eq!(fpmath::sinh(-x), -actual_sin1);
        assert_total_eq!(fpmath::cosh(-x), actual_cos1);
        assert_total_eq!(fpmath::sinh_cosh(-x), (-actual_sin2, actual_cos2));

        check_result(x, actual_sin1, expected_sin, 0.51, &mut max_sin1_error);
        check_result(x, actual_sin2, expected_sin, 0.51, &mut max_sin2_error);
        check_result(x, actual_cos1, expected_cos, 0.51, &mut max_cos1_error);
        check_result(x, actual_cos2, expected_cos, 0.51, &mut max_cos2_error);
    });
    eprintln!("max sinh1 error = {max_sin1_error}");
    eprintln!("max sinh2 error = {max_sin2_error}");
    eprintln!("max cosh1 error = {max_cos1_error}");
    eprintln!("max cosh2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_tanh() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tanh(f64::from(x));
        let actual = fpmath::tanh(x);
        assert_total_eq!(fpmath::tanh(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=9 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..10000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for arg in 1..=1000 {
        f(arg as f32);
    }
}
