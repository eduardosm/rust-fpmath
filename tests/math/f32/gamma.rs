use super::{calc_error_ulp, mk_normal, mk_subnormal, purify};
use crate::create_prng;

#[test]
fn test_gamma() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::gamma(f64::from(x));
        let actual = fpmath::gamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "gamma({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_ln_gamma() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let (expected, expected_sign) = fpmath::ln_gamma(f64::from(x));
        let (actual, actual_sign) = fpmath::ln_gamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert_eq!(expected_sign, actual_sign);
        assert!(
            err < 0.51,
            "ln_gamma({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error >= 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    // Exhaustive test of all subnormal numbers
    for m in 0..(1 << 23) {
        f(mk_subnormal(m, false));
        f(mk_subnormal(m, true));
    }

    for e in -126..=127 {
        for _ in 0..6000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
            f(mk_normal(m, e, true));
        }
    }

    // Problematic range
    for e in 1..=2 {
        for m in 0..(1 << 23) {
            f(mk_normal(m, e, true));
        }
    }

    for i in 0..20000 {
        let x = purify((i as f32) / 100.0);
        f(x);
        f(-x);
    }

    // Some roots of ln_gamma (i.e., where abs(gamma(x)) = 1 and ln_gamma(x) = 0)
    let roots = [
        1.0, 2.0, -2.4570246, -2.7476826, -3.143581, -3.9552946, -4.039362, -4.9915447, -5.0082183,
        -5.9986076,
    ];

    for root in roots {
        let root = root as f32;
        f(root);

        for bump in 1..=2000 {
            let x = f32::from_bits(root.to_bits() + bump);
            f(x);

            let x = f32::from_bits(root.to_bits() - bump);
            f(x);
        }
    }
}
