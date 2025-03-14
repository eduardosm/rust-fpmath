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

fn mkfloat(m: u64, e: i16, s: bool) -> f64 {
    let m = m >> (64 - 52);
    let e = u64::from((e + 1023) as u16) << 52;
    let s = u64::from(s) << 63;
    f64::from_bits(m | e | s)
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

    match expected.classify() {
        std::num::FpCategory::Nan => {
            if actual.is_nan() {
                0.0
            } else {
                f64::INFINITY
            }
        }
        std::num::FpCategory::Infinite => {
            if actual.is_infinite() && actual.is_sign_positive() == expected.is_sign_positive() {
                0.0
            } else {
                f64::INFINITY
            }
        }
        std::num::FpCategory::Subnormal => unreachable!(),
        _ if actual.is_infinite() => {
            if expected.get_exp().is_some_and(|e| e > 1023)
                && actual.is_sign_positive() == expected.is_sign_positive()
            {
                0.0
            } else {
                f64::INFINITY
            }
        }
        _ => {
            let exp = expected
                .get_exp()
                .map(|e| (e - 1).max(-1022))
                .unwrap_or(-1022);
            let dif = (expected - actual).abs() >> (exp - 52);
            dif.to_f64()
        }
    }
}

// Workaround X87 compiler bugs
fn purify(x: f64) -> f64 {
    std::hint::black_box(x)
}

fn purify2((x, y): (f64, f64)) -> (f64, f64) {
    (std::hint::black_box(x), std::hint::black_box(y))
}
