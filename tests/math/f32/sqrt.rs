use super::{mk_normal, mk_subnormal};
use crate::create_prng;

#[test]
fn test_sqrt() {
    test_with(|x| {
        let actual = fpmath::sqrt(x);
        let expected = rug::Float::with_val(24, x).sqrt();
        assert_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    // Exhaustive test of all mantissas
    for e in 0..=1 {
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
