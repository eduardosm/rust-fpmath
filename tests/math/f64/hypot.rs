use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, select_threshold, RUG_PREC};
use crate::create_prng;

#[test]
fn test_hypot() {
    let mut max_error: f64 = 0.0;
    test_with(|x, y| {
        let expected = rug::Float::with_val(RUG_PREC, x).hypot(&rug::Float::with_val(RUG_PREC, y));
        let actual = fpmath::hypot(x, y);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "hypot({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max hypot error = {max_error}");
    assert!(max_error > 0.49);
}

fn test_with(mut f: impl FnMut(f64, f64)) {
    let mut rng = create_prng();

    for ex in -1022..=1023 {
        if matches!(ex, -900..=900) && (ex & 3) != 3 {
            continue; // speed up tests
        }
        for ey in -1022..=1023 {
            for _ in 0..5 {
                let mx = rng.random::<u64>();
                let sx = rng.random::<bool>();
                let my = rng.random::<u64>();
                let sy = rng.random::<bool>();
                f(mkfloat(mx, ex, sx), mkfloat(my, ey, sy));
            }
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u64>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, e, sx), mkfloat(my, e, sy));

            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u64>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, 0, sx), mkfloat(my, e, sy));

            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u64>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, e, sx), mkfloat(my, 0, sy));
        }
    }
}
