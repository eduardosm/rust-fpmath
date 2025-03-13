use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, purify};
use crate::data::create_prng;

#[test]
fn test_tgamma() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tgamma(f64::from(x));
        let actual = fpmath::tgamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = if x < 0.5 { 1.9 } else { 0.9 };
        assert!(
            err < threshold,
            "tgamma({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max tgamma error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_lgamma() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let (expected, expected_sign) = fpmath::lgamma(f64::from(x));
        let (actual, actual_sign) = fpmath::lgamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = if (-5.0..=-2.0).contains(&x) {
            // FIXME
            200.0
        } else if (0.5..=7.0).contains(&x) {
            1.5
        } else {
            1.9
        };
        assert_eq!(expected_sign, actual_sign);
        assert!(
            err < threshold,
            "lgamma({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max lgamma error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        for _ in 0..6000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
            f(mkfloat(m, e, true));
        }
    }

    for i in 0..20000 {
        let x = purify((i as f32) / 100.0);
        f(x);
        f(-x);
    }
}
