use super::{calc_error_ulp, mk_normal, mk_subnormal};
use crate::create_prng;

#[test]
fn test_cbrt() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::cbrt(f64::from(x));
        let actual = fpmath::cbrt(x);
        assert_total_eq!(fpmath::cbrt(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "cbrt({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.4999);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    // Exhaustive test of all mantissas
    for e in 0..=2 {
        for m in 0..(1 << 23) {
            f(mk_normal(m, e, false));
        }
    }

    // Exhaustive test of all subnormal numbers
    for m in 0..(1 << 23) {
        f(mk_subnormal(m, false));
    }

    for e in -126..=127 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..1000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for arg in 1..=500_000 {
        f(arg as f32);
    }

    f(f32::MIN_POSITIVE);
    f(f32::MAX);
}
