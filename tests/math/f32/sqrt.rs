use rand::Rng as _;

use super::mkfloat;
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

    for e in -126..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(u32::MAX, e, false));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f32);
    }

    f(f32::MIN_POSITIVE);
    f(f32::MAX);

    // subnormals
    for i in 0..23 {
        f(f32::from_bits(1 << i));
        f(f32::from_bits((1 << (i + 1)) - 1));
    }
}
