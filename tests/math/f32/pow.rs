use rand::RngExt as _;

use super::{calc_error_ulp, mk_normal, mk_subnormal, purify};
use crate::create_prng;

#[test]
fn test_pow() {
    let mut max_error: f32 = 0.0;
    test_pow_with(|x, y| {
        let expected = fpmath::pow(f64::from(x), f64::from(y));
        let actual = fpmath::pow(x, y);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "pow({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_pow_with(mut f: impl FnMut(f32, f32)) {
    let mut rng = create_prng();

    // x = sx * mx * 2^ex
    // log2(|x|) = log2(mx) + ex
    // ex <= log2(|x|) <= ex + 1

    // MIN <= |x|^y <= MAX
    // log2(MIN) / log2(|x|) <= y <= log2(MAX) / log2(|x|)

    for ex in -126..=127 {
        let (min_y, max_y) = if ex == 0 {
            (-126, 127)
        } else {
            let a = 127 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-126), (a.max(b) + 3).min(127))
        };

        for yi in min_y..=max_y {
            for _ in 0..100 {
                let mx = super::gen_mantissa(&mut rng);
                let sx = false;
                let x = mk_normal(mx, ex, sx);

                let y = (rng.random::<f32>() - 0.5) + (yi as f32);
                f(x, purify(y));
            }
        }
    }

    for ex in -22..=-1 {
        for ey in 1..=127 {
            for _ in 0..2000 {
                let mx = super::gen_mantissa(&mut rng);
                let my = super::gen_mantissa(&mut rng);

                f(
                    purify(1.0 + mk_normal(mx, ex, false)),
                    mk_normal(my, ey, false),
                );
                f(
                    purify(1.0 + mk_normal(mx, ex, false)),
                    mk_normal(my, ey, true),
                );
                f(
                    purify(1.0 + mk_normal(mx, ex, true)),
                    mk_normal(my, ey, false),
                );
                f(
                    purify(1.0 + mk_normal(mx, ex, true)),
                    mk_normal(my, ey, true),
                );
            }
        }
    }

    for _ in 0..100_000 {
        let mx = super::gen_mantissa(&mut rng);
        let x = mk_subnormal(mx, false);
        let my = super::gen_mantissa(&mut rng);
        let y = mk_normal(my, -1, false);

        f(x, y);
        f(x, -y);
    }
}

#[test]
fn test_powi() {
    let mut max_error: f32 = 0.0;
    test_powi_with(|x, y| {
        let expected = fpmath::powi(f64::from(x), y);
        let actual = fpmath::powi(x, y);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(
            err < 0.51,
            "powi({x:e}, {y}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_powi_with(mut f: impl FnMut(f32, i32)) {
    let mut rng = create_prng();

    for ex in -126..=127 {
        let (min_y, max_y) = if ex == 0 {
            (-126, 127)
        } else {
            let a = 127 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-126), (a.max(b) + 3).min(127))
        };

        for y in min_y..=max_y {
            for _ in 0..100 {
                let mx = super::gen_mantissa(&mut rng);
                let x = mk_normal(mx, ex, false);

                f(x, i32::from(y));
            }
        }
    }

    for ex in -22..=-1 {
        for i in (1..=31).rev() {
            for _ in 0..2000 {
                let mx = super::gen_mantissa(&mut rng);
                let y = ((rng.random::<u32>() | 0x8000_0000) >> i) as i32;

                let xp = purify(1.0 + mk_normal(mx, ex, false));
                let xn = purify(1.0 + mk_normal(mx, ex, true));

                f(xp, y);
                f(xp, -y);
                f(xn, y);
                f(xn, -y);
            }
        }
    }

    for _ in 0..100000 {
        let mx = super::gen_mantissa(&mut rng);
        let x = mk_subnormal(mx, false);

        f(x, 1);
        f(x, -1);
    }
}
