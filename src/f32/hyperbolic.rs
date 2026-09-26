use crate::generic::{round_fi, scalbn};
use crate::traits::Float as _;

// GENERATE: consts f32 LOG2_E
const LOG2_E: f32 = f32::from_bits(0x3FB8AA3B); // 1.442695e0

// GENERATE: consts f64 LN_2
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1

impl crate::generic::Hyperbolic for f32 {
    fn sinh_finite(x: Self) -> Self {
        if x.exponent() >= 7 {
            f32::INFINITY.copysign(x)
        } else {
            let (sinh, _) = sinh_cosh_inner(x);
            sinh as f32
        }
    }

    fn cosh_finite(x: Self) -> Self {
        if x.exponent() >= 7 {
            f32::INFINITY
        } else {
            let (_, cosh) = sinh_cosh_inner(x);
            cosh as f32
        }
    }

    fn sinh_cosh_finite(x: Self) -> (Self, Self) {
        if x.exponent() >= 7 {
            (f32::INFINITY.copysign(x), f32::INFINITY)
        } else {
            let (sinh, cosh) = sinh_cosh_inner(x);
            (sinh as f32, cosh as f32)
        }
    }

    fn tanh_finite(x: Self) -> Self {
        if x.exponent() >= 5 {
            1.0f32.copysign(x)
        } else {
            let (sinh, cosh) = sinh_cosh_inner(x);
            (sinh / cosh) as f32
        }
    }
}

fn sinh_cosh_inner(x: f32) -> (f64, f64) {
    #[inline]
    fn exp_m1_core(r: f64) -> (f64, f64) {
        // exp(r) = exp(r / s)^s
        // exp(-r) = exp(-r / s)^s
        // s = 2^4

        let rs = r * (1.0 / 16.0);

        // GENERATE: exp_m1_poly f64 3 -0.0216608 0.0216608
        const K2: f64 = f64::from_bits(0x3FDFFFFFFFDD2FFD); // 4.999999998733527e-1
        const K3: f64 = f64::from_bits(0x3FC55570814769D0); // 1.6666990578126084e-1
        const K4: f64 = f64::from_bits(0x3FA5557441B2698B); // 4.166758825268165e-2

        let rs2 = rs * rs;
        let mrs = -rs;
        let mut t1 = rs + horner!(rs2, rs, [K2, K3, K4]);
        let mut t2 = mrs + horner!(rs2, mrs, [K2, K3, K4]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..4 {
            t1 = t1 * t1 + 2.0 * t1;
            t2 = t2 * t2 + 2.0 * t2;
        }

        (t1, t2)
    }

    if x.exponent() >= 7 {
        return (f64::INFINITY.set_sign(x.is_sign_negative()), f64::INFINITY);
    }

    // x = k*ln(2) + r
    // k is an integer
    // |r| <= 0.5*ln(2)
    let (kf, k) = round_fi(x * LOG2_E);
    let r = f64::from(x) - f64::from(kf) * LN_2;

    // t1 = exp(r) - 1
    // t2 = exp(-r) - 1
    let (t1a, t1b) = exp_m1_core(r);

    // sinh(x) = (exp(x) - exp(-x)) / 2
    // cosh(x) = (exp(x) + exp(-x)) / 2

    if k == 0 {
        let s = 0.5 * (t1a - t1b);
        let c = 0.5 * (t1a + t1b) + 1.0;
        (s, c)
    } else {
        let t2a = scalbn(t1a + 1.0, k - 1);
        let t2b = scalbn(t1b + 1.0, -k - 1);

        let s = t2a - t2b;
        let c = t2a + t2b;
        (s, c)
    }
}
