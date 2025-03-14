use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, select_threshold, RUG_PREC};
use crate::create_prng;

#[test]
fn test_pow() {
    let mut max_error: f64 = 0.0;
    test_pow_with(|x, y| {
        let bigx = rug::Float::with_val(53, x);
        let bigy = rug::Float::with_val(53, y);
        let expected = rug::Float::with_val(RUG_PREC, rug::ops::Pow::pow(&bigx, &bigy));
        let actual = fpmath::pow(x, y);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "pow({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max pow error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_pow_with(mut f: impl FnMut(f64, f64)) {
    let mut rng = create_prng();

    // x = sx * mx * 2^ex
    // log2(|x|) = log2(mx) + ex
    // ex <= log2(|x|) <= ex + 1

    // MIN <= |x|^y <= MAX
    // log2(MIN) / log2(|x|) <= y <= log2(MAX) / log2(|x|)

    for ex in -1022..=1023 {
        let (min_y, max_y) = if ex == 0 {
            (-1022, 1023)
        } else {
            let a = 1023 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-1022), (a.max(b) + 3).min(1023))
        };

        for yi in min_y..=max_y {
            for _ in 0..50 {
                let mx = rng.random::<u64>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                let y = (rng.random::<f64>() - 0.5) + (yi as f64);
                f(x, y);
            }
        }
    }

    for ex in -51..=-1 {
        for ey in 1..=1023 {
            for _ in 0..100 {
                let mx = rng.random::<u64>();
                let sx = rng.random::<bool>();
                let my = rng.random::<u64>();
                let sy = rng.random::<bool>();
                f(1.0 + mkfloat(mx, ex, sx), mkfloat(my, ey, sy));
            }
        }
    }
}

#[test]
fn test_powi() {
    let mut max_error: f64 = 0.0;
    test_powi_with(|x, y| {
        let bigx = rug::Float::with_val(53, x);
        let expected = rug::Float::with_val(RUG_PREC, rug::ops::Pow::pow(&bigx, y));
        let actual = fpmath::powi(x, y);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "powi({x:e}, {y}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max pow error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_powi_with(mut f: impl FnMut(f64, i32)) {
    let mut rng = create_prng();

    for ex in -1022..=1023 {
        let (min_y, max_y) = if ex == 0 {
            (-1022, 1023)
        } else {
            let a = 1023 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-1022), (a.max(b) + 3).min(1023))
        };

        for y in min_y..=max_y {
            for _ in 0..100 {
                let mx = rng.random::<u64>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                f(x, i32::from(y));
            }
        }
    }

    for ex in -51..=-1 {
        for i in (1..=31).rev() {
            for _ in 0..1000 {
                let mx = rng.random::<u64>();
                let sx = rng.random::<bool>();
                let x = 1.0 + mkfloat(mx, ex, sx);
                let y = ((rng.random::<u32>() | 0x8000_0000) >> i) as i32;

                f(x, y);
                f(x, -y);
            }
        }
    }
}
