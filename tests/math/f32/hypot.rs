use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, select_threshold};
use crate::create_prng;

#[test]
fn test_hypot() {
    let mut max_error: f32 = 0.0;
    test_with(|x, y| {
        let expected = fpmath::hypot(f64::from(x), f64::from(y));
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

fn test_with(mut f: impl FnMut(f32, f32)) {
    let mut rng = create_prng();

    for ex in -126..=127 {
        for ey in -126..=127 {
            for _ in 0..5 {
                let mx = rng.random::<u32>();
                let sx = rng.random::<bool>();
                let my = rng.random::<u32>();
                let sy = rng.random::<bool>();
                f(mkfloat(mx, ex, sx), mkfloat(my, ey, sy));
            }
        }
    }

    for e in -126..=127 {
        for _ in 0..5000 {
            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u32>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, e, sx), mkfloat(my, e, sy));

            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u32>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, 0, sx), mkfloat(my, e, sy));

            let mx = rng.random::<u32>();
            let sx = rng.random::<bool>();
            let my = rng.random::<u32>();
            let sy = rng.random::<bool>();
            f(mkfloat(mx, e, sx), mkfloat(my, 0, sy));
        }
    }
}
