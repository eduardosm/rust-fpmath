use super::mkfloat;
use crate::create_prng;

#[test]
fn test_round() {
    test_round_with(|arg| {
        let expected = rug::Float::with_val(128, arg).round();
        let actual = fpmath::round(arg);

        assert!(expected == actual, "round({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_floor() {
    test_round_with(|arg| {
        let expected = rug::Float::with_val(128, arg).floor();
        let actual = fpmath::floor(arg);

        assert!(expected == actual, "floor({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_ceil() {
    test_round_with(|arg| {
        let expected = rug::Float::with_val(128, arg).ceil();
        let actual = fpmath::ceil(arg);

        assert!(expected == actual, "ceil({arg:e}) = {actual:e}");
    });
}

#[test]
fn test_trunc() {
    test_round_with(|arg| {
        let expected = rug::Float::with_val(128, arg).trunc();
        let actual = fpmath::trunc(arg);

        assert!(expected == actual, "trunc({arg:e}) = {actual:e}");
    });
}

fn test_round_with(test_f: fn(f64)) {
    use rand::Rng as _;

    let mut rng = create_prng();

    for e in -1022..=1023 {
        test_f(mkfloat(0, e, false));
        test_f(mkfloat(0, e, true));
        test_f(mkfloat(u64::MAX, e, false));
        test_f(mkfloat(u64::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.random::<u64>();
            test_f(mkfloat(m, e, true));
            test_f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=100_000 {
        let arg = arg as f64;
        test_f(arg);
        test_f(-arg);
        test_f(arg + 0.25);
        test_f(-arg + 0.25);
        test_f(arg + 0.5);
        test_f(-arg + 0.5);
        test_f(arg + 0.75);
        test_f(-arg + 0.75);
    }

    for e in 0..=52 {
        for delta in -1000..=1000 {
            let arg = mkfloat(0, e, false) + delta as f64;
            test_f(arg);
            test_f(-arg);
            test_f(arg + 0.25);
            test_f(-arg + 0.25);
            test_f(arg + 0.5);
            test_f(-arg + 0.5);
            test_f(arg + 0.75);
            test_f(-arg + 0.75);
        }
    }
}
