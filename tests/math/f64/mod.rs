mod cbrt;
mod exp;
mod gamma;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod round;
mod sqrt;
mod trigonometric;

const MAX_MANTISSA: u64 = (1 << 52) - 1;

fn gen_mantissa(rng: &mut impl rand::RngExt) -> u64 {
    rng.random::<u64>() >> (64 - 52)
}

fn mk_normal(m: u64, e: i16, s: bool) -> f64 {
    assert!(m < (1 << 52));
    assert!(matches!(e, -1022..=1023));
    let e = u64::from((e + 1023) as u16) << 52;
    let s = u64::from(s) << 63;
    f64::from_bits(m | e | s)
}

fn mk_subnormal(m: u64, s: bool) -> f64 {
    assert!(m < (1 << 52));
    let s = u64::from(s) << 63;
    f64::from_bits(m | s)
}

fn select_threshold(actual: f64, normal_th: f64, subnormal_th: f64) -> f64 {
    if actual == 0.0 || actual.is_subnormal() {
        subnormal_th
    } else {
        normal_th
    }
}

const RUG_PREC: u32 = 53 + 20;

fn calc_error_ulp(actual: f64, expected: rug::Float) -> f64 {
    let actual = purify(actual);

    if expected.is_nan() {
        if actual.is_nan() { 0.0 } else { f64::INFINITY }
    } else if expected > f64::MAX {
        if actual == f64::INFINITY {
            0.0
        } else {
            f64::INFINITY
        }
    } else if expected < f64::MIN {
        if actual == f64::NEG_INFINITY {
            0.0
        } else {
            f64::INFINITY
        }
    } else if actual.is_infinite() {
        f64::INFINITY
    } else {
        let exp = expected
            .get_exp()
            .map(|e| (e - 1).max(-1022))
            .unwrap_or(-1022);
        let dif = (expected - actual).abs() >> (exp - 52);
        dif.to_f64()
    }
}

// Workaround X87 compiler bugs
fn purify(x: f64) -> f64 {
    std::hint::black_box(x)
}

impl crate::TotalEq for f64 {
    fn total_eq(&self, other: &Self) -> bool {
        self.to_bits() == other.to_bits()
    }
}
