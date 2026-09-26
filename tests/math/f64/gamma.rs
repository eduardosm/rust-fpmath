use super::{RUG_PREC, calc_error_ulp, mk_normal, mk_subnormal, purify};
use crate::create_prng;

#[test]
fn test_gamma() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).gamma();
        let actual = fpmath::gamma(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "gamma({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error >= 0.5);
}

#[test]
fn test_ln_gamma() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let (expected, ord) = rug::Float::with_val(RUG_PREC, x).ln_abs_gamma();
        let expected_sign = if x < 0.0 && x.fract() == 0.0 {
            0
        } else {
            ord as i8
        };
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
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        for _ in 0..3000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
            f(mk_normal(m, e, true));
        }
    }

    for _ in 0..100000 {
        let m = super::gen_mantissa(&mut rng);
        f(mk_subnormal(m, false));
        f(mk_subnormal(m, true));
    }

    for i in 0..51 {
        let m = 1 << i;
        f(mk_subnormal(m, false));
        f(mk_subnormal(m, true));
    }

    // Problematic range
    for e in 1..=2 {
        for _ in 0..2_000_000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, true));
        }
    }

    for i in 0..20000 {
        let x = purify((i as f64) / 100.0);
        f(x);
        f(-x);
    }

    // Roots of ln_gamma (i.e., where abs(gamma(x)) = 1 and ln_gamma(x) = 0)
    let roots = [
        1.0,
        2.0,
        -2.4570247382208006,
        -2.7476826467274127,
        -3.14358088834998,
        -3.955294284858598,
        -4.039361839740537,
        -4.991544640560048,
        -5.0082181683225935,
        -5.998607480080875,
        -6.001385294453155,
        -6.999801507890638,
        -7.000198333407325,
        -7.999975197095821,
        -8.000024800270682,
        -8.999997244250977,
        -9.000002755714823,
        -9.99999972442663,
        -10.000000275573013,
        -10.99999997494789,
        -11.000000025052106,
        -11.999999997912324,
        -12.000000002087676,
        -12.99999999983941,
        -13.00000000016059,
        -13.99999999998853,
        -14.00000000001147,
        -14.999999999999236,
        -15.000000000000764,
        -15.999999999999952,
        -16.000000000000046,
        -16.999999999999996,
        -17.000000000000004,
    ];

    for root in roots {
        f(root);

        for bump in 1..=10_000 {
            let x = f64::from_bits(root.to_bits() + bump);
            f(x);

            let x = f64::from_bits(root.to_bits() - bump);
            f(x);
        }
    }
}
