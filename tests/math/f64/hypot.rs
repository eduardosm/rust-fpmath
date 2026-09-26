use super::{RUG_PREC, calc_error_ulp, mk_normal, mk_subnormal};
use crate::create_prng;

#[test]
fn test_hypot() {
    let mut max_error: f64 = 0.0;
    test_with(|x, y| {
        let expected = rug::Float::with_val(RUG_PREC, x).hypot(&rug::Float::with_val(RUG_PREC, y));
        let actual = fpmath::hypot(x, y);
        assert_total_eq!(fpmath::hypot(-x, y), actual);
        assert_total_eq!(fpmath::hypot(x, -y), actual);
        assert_total_eq!(fpmath::hypot(-x, -y), actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "hypot({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.4999);
}

fn test_with(mut f: impl FnMut(f64, f64)) {
    let mut rng = create_prng();

    for ex in -1022..=1023 {
        if matches!(ex, -900..=900) && (ex & 3) != 3 {
            continue; // speed up tests
        }
        for ey in -1022..=1023 {
            for _ in 0..5 {
                let mx = super::gen_mantissa(&mut rng);
                let my = super::gen_mantissa(&mut rng);
                f(mk_normal(mx, ex, false), mk_normal(my, ey, false));
            }
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let mx = super::gen_mantissa(&mut rng);
            let my = super::gen_mantissa(&mut rng);
            f(mk_normal(mx, e, false), mk_normal(my, e, false));

            let mx = super::gen_mantissa(&mut rng);
            let my = super::gen_mantissa(&mut rng);
            f(mk_normal(mx, 0, false), mk_normal(my, e, false));

            let mx = super::gen_mantissa(&mut rng);
            let my = super::gen_mantissa(&mut rng);
            f(mk_normal(mx, e, false), mk_normal(my, 0, false));
        }
    }

    for ix in 0..52 {
        let x = mk_subnormal(1 << ix, false);

        f(x, 0.0);
        f(0.0, x);

        for iy in 0..52 {
            let y = mk_subnormal(1 << iy, false);
            f(x, y);
        }
    }

    for _ in 0..100_000 {
        let mx = super::gen_mantissa(&mut rng);
        let x = mk_subnormal(mx, false);
        let my = super::gen_mantissa(&mut rng);
        let y = mk_subnormal(my, false);

        f(x, y);
        f(x, 0.0);
        f(0.0, y);
    }
}
