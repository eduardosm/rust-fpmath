use super::mk_normal;
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
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..10000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
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
