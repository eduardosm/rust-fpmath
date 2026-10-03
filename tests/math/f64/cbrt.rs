use super::{RUG_PREC, check_result, mk_normal, mk_subnormal};
use crate::create_prng;

#[test]
fn test_cbrt() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).cbrt();
        let actual = fpmath::cbrt(x);
        assert_total_eq!(fpmath::cbrt(-x), -actual);

        check_result(x, actual, expected, 0.51, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error >= 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..5000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..5000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f64);
    }

    f(f64::MIN_POSITIVE);
    f(f64::MAX);

    // subnormals
    for i in 0..52 {
        f(mk_subnormal(1 << i, false));
        f(mk_subnormal((1 << (i + 1)) - 1, false));
    }
}
