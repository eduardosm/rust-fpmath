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

const MAX_MANTISSA: u32 = (1 << 23) - 1;

fn gen_mantissa(rng: &mut impl rand::RngExt) -> u32 {
    rng.random::<u32>() >> (32 - 23)
}

fn mk_normal(m: u32, e: i16, s: bool) -> f32 {
    assert!(matches!(e, -126..=127));
    let e = u32::from((e + 127) as u16) << 23;
    let s = u32::from(s) << 31;
    f32::from_bits(m | e | s)
}

fn mk_subnormal(m: u32, s: bool) -> f32 {
    assert!(m < (1 << 23));
    let s = u32::from(s) << 31;
    f32::from_bits(m | s)
}

fn calc_error_ulp(actual: f32, expected: f64) -> f32 {
    let actual = purify(actual);

    if expected.is_nan() {
        if actual.is_nan() { 0.0 } else { f32::INFINITY }
    } else if actual.is_nan() {
        f32::INFINITY
    } else if expected > f64::from(f32::MAX) {
        if actual == f32::INFINITY {
            0.0
        } else {
            f32::INFINITY
        }
    } else if expected < f64::from(f32::MIN) {
        if actual == f32::NEG_INFINITY {
            0.0
        } else {
            f32::INFINITY
        }
    } else if actual.is_infinite() {
        f32::INFINITY
    } else {
        let exp = if expected == 0.0 {
            -126
        } else {
            (fpmath::frexp(expected).1 - 1).max(-126)
        };
        let dif = fpmath::scalbn((expected - f64::from(actual)).abs(), 23 - exp);
        dif as f32
    }
}

// Workaround X87 compiler bugs
fn purify(x: f32) -> f32 {
    std::hint::black_box(x)
}

impl crate::TotalEq for f32 {
    fn total_eq(&self, other: &Self) -> bool {
        self.to_bits() == other.to_bits()
    }
}
