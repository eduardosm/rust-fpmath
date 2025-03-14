use rand::Rng as _;

use super::mkfloat;
use crate::create_prng;

#[test]
fn test_sqrt() {
    test_with(|x| {
        let actual = fpmath::sqrt(x);
        let expected = rug::Float::with_val(53, x).sqrt();
        assert_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -200..=200 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..5000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f64);
    }

    f(f64::MIN_POSITIVE);
    f(f64::MAX);

    // subnormals
    for i in 0..52 {
        f(f64::from_bits(1 << i));
        f(f64::from_bits((1 << (i + 1)) - 1));
    }
}
