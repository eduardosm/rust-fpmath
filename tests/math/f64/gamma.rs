use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, RUG_PREC};
use crate::create_prng;

#[test]
fn test_tgamma() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).gamma();
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
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let (expected, ord) = rug::Float::with_val(RUG_PREC, x).ln_abs_gamma();
        let expected_sign = if x < 0.0 && x.fract() == 0.0 {
            0
        } else {
            ord as i8
        };
        let (actual, actual_sign) = fpmath::lgamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = if (-5.0..=-2.0).contains(&x) {
            // FIXME
            50.0
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

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        for _ in 0..3000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
            f(mkfloat(m, e, true));
        }
    }

    for i in 0..20000 {
        let x = (i as f64) / 100.0;
        f(x);
        f(-x);
    }
}
