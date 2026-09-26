use crate::generic::{reduce_90_deg, reduce_half_revs, reduce_pi_2_large, round_fi, scalbn};
use crate::traits::Float as _;

// GENERATE: consts f64 FRAC_PI_180 PI
const FRAC_PI_180: f64 = f64::from_bits(0x3F91DF46A2529D39); // 1.7453292519943295e-2
const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0

impl crate::generic::Trigonometric for f32 {
    fn sin_finite(x: Self) -> Self {
        let (n, y) = reduce_pi_2(x);
        match n {
            0 => sin_core(y),
            1 => cos_core(y),
            2 => -sin_core(y),
            3 => -cos_core(y),
            _ => unreachable!(),
        }
    }

    fn cos_finite(x: Self) -> Self {
        let (n, y) = reduce_pi_2(x);
        match n {
            0 => cos_core(y),
            1 => -sin_core(y),
            2 => -cos_core(y),
            3 => sin_core(y),
            _ => unreachable!(),
        }
    }

    fn sin_cos_finite(x: Self) -> (Self, Self) {
        let (n, y) = reduce_pi_2(x);
        let sin = sin_core(y);
        let cos = cos_core(y);
        let (sin, cos) = match n {
            0 => (sin, cos),
            1 => (cos, -sin),
            2 => (-sin, -cos),
            3 => (-cos, sin),
            _ => unreachable!(),
        };
        (sin, cos)
    }

    fn tan_finite(x: Self) -> Self {
        let (n, y) = reduce_pi_2(x);
        tan_core(y, (n & 1) != 0)
    }

    fn sind_finite(x: Self) -> Self {
        if x.raw_exp() == 0 {
            (f64::from(x) * FRAC_PI_180) as f32
        } else {
            let (n, y) = reduce_90_deg(x);
            if y == 0.0 {
                match n {
                    0 => 0.0f32.copysign(x),
                    1 => 1.0,
                    2 => 0.0f32.copysign(x),
                    3 => -1.0,
                    _ => unreachable!(),
                }
            } else {
                let z = f64::from(y) * FRAC_PI_180;
                match n {
                    0 => sin_core(z),
                    1 => cos_core(z),
                    2 => -sin_core(z),
                    3 => -cos_core(z),
                    _ => unreachable!(),
                }
            }
        }
    }

    fn cosd_finite(x: Self) -> Self {
        let (n, y) = reduce_90_deg(x);
        if y == 0.0 {
            match n {
                0 => 1.0,
                1 => 0.0,
                2 => -1.0,
                3 => 0.0,
                _ => unreachable!(),
            }
        } else {
            let z = f64::from(y) * FRAC_PI_180;
            match n {
                0 => cos_core(z),
                1 => -sin_core(z),
                2 => -cos_core(z),
                3 => sin_core(z),
                _ => unreachable!(),
            }
        }
    }

    fn sind_cosd_finite(x: Self) -> (Self, Self) {
        if x.raw_exp() == 0 {
            ((f64::from(x) * FRAC_PI_180) as f32, 1.0)
        } else {
            let (n, y) = reduce_90_deg(x);
            if y == 0.0 {
                match n {
                    0 => (0.0f32.copysign(x), 1.0),
                    1 => (1.0, 0.0),
                    2 => (0.0f32.copysign(x), -1.0),
                    3 => (-1.0, 0.0),
                    _ => unreachable!(),
                }
            } else {
                let z = f64::from(y) * FRAC_PI_180;
                let sin = sin_core(z);
                let cos = cos_core(z);
                match n {
                    0 => (sin, cos),
                    1 => (cos, -sin),
                    2 => (-sin, -cos),
                    3 => (-cos, sin),
                    _ => unreachable!(),
                }
            }
        }
    }

    fn tand_finite(x: Self) -> Self {
        if x.raw_exp() == 0 {
            (f64::from(x) * FRAC_PI_180) as f32
        } else {
            let (n, y) = reduce_90_deg(x);
            let inv = (n & 1) != 0;
            if y == 0.0 {
                if inv {
                    f32::INFINITY.set_sign(n == 3)
                } else {
                    0.0f32.set_sign(x.is_sign_negative() ^ (n == 2))
                }
            } else {
                let z = f64::from(y) * FRAC_PI_180;
                tan_core(z, inv)
            }
        }
    }

    fn sinpi_finite(x: Self) -> Self {
        if x.raw_exp() == 0 {
            (f64::from(x) * PI) as f32
        } else {
            let (n, y) = reduce_half_revs(x);
            if y == 0.0 {
                match n {
                    0 => 0.0f32.copysign(x),
                    1 => 1.0,
                    2 => 0.0f32.copysign(x),
                    3 => -1.0,
                    _ => unreachable!(),
                }
            } else {
                let z = f64::from(y) * PI;
                match n {
                    0 => sin_core(z),
                    1 => cos_core(z),
                    2 => -sin_core(z),
                    3 => -cos_core(z),
                    _ => unreachable!(),
                }
            }
        }
    }

    fn cospi_finite(x: Self) -> Self {
        let (n, y) = reduce_half_revs(x);
        if y == 0.0 {
            match n {
                0 => 1.0,
                1 => 0.0,
                2 => -1.0,
                3 => 0.0,
                _ => unreachable!(),
            }
        } else {
            let z = f64::from(y) * PI;
            match n {
                0 => cos_core(z),
                1 => -sin_core(z),
                2 => -cos_core(z),
                3 => sin_core(z),
                _ => unreachable!(),
            }
        }
    }

    fn sinpi_cospi_finite(x: Self) -> (Self, Self) {
        if x.raw_exp() == 0 {
            ((f64::from(x) * PI) as f32, 1.0)
        } else {
            let (n, y) = reduce_half_revs(x);
            if y == 0.0 {
                match n {
                    0 => (0.0f32.copysign(x), 1.0),
                    1 => (1.0, 0.0),
                    2 => (0.0f32.copysign(x), -1.0),
                    3 => (-1.0, 0.0),
                    _ => unreachable!(),
                }
            } else {
                let z = f64::from(y) * PI;
                let sin = sin_core(z);
                let cos = cos_core(z);
                match n {
                    0 => (sin, cos),
                    1 => (cos, -sin),
                    2 => (-sin, -cos),
                    3 => (-cos, sin),
                    _ => unreachable!(),
                }
            }
        }
    }

    fn tanpi_finite(x: Self) -> Self {
        if x.raw_exp() == 0 {
            (f64::from(x) * PI) as f32
        } else {
            let (n, y) = reduce_half_revs(x);
            let inv = (n & 1) != 0;
            if y == 0.0 {
                if inv {
                    f32::INFINITY.set_sign(n == 3)
                } else {
                    0.0f32.set_sign(x.is_sign_negative() ^ (n == 2))
                }
            } else {
                let z = f64::from(y) * PI;
                tan_core(z, inv)
            }
        }
    }
}

/// Reduces the angle argument `x`, returning `(n, y)` such that:
/// * `|y| <= π/4`
/// * `0 <= n <= 3`
/// * `x = 2*π*M + π/2*n + y`
/// * `M` is an integer
fn reduce_pi_2(x: f32) -> (u8, f64) {
    let x_orig = x;

    // GENERATE: consts f32 FRAC_PI_4 FRAC_2_PI
    const FRAC_PI_4: f32 = f32::from_bits(0x3F490FDB); // 7.853982e-1
    const FRAC_2_PI: f32 = f32::from_bits(0x3F22F983); // 6.3661975e-1

    // GENERATE: reduce_pi_2::medium_consts f32 f64
    const FRAC_PI_2_P0: f32 = f32::from_bits(0x3FC90FDA); // 1.5707963e0
    const FRAC_PI_2_P0EX: f64 = f64::from_bits(0x3E74442D18469899); // 7.549789954891882e-8

    if x.abs() <= FRAC_PI_4 {
        // reduction not needed
        return (0, x.into());
    }

    let p = x * FRAC_2_PI;
    if p.exponent() <= 18 {
        let (f_n, n) = round_fi(p);

        let x = f64::from(x);
        let f_n = f64::from(f_n);

        let r = x - f_n * f64::from(FRAC_PI_2_P0);
        let y = r - f_n * FRAC_PI_2_P0EX;

        (n as u8 & 3, y)
    } else {
        let mant = x_orig.mant();
        let x_chunks = [mant];
        let e0 = x_orig.exponent() - 23;
        let jk = 4;

        let mut qp: [u64; 20] = [0; 20];
        let (ih, jz, n, qe) = reduce_pi_2_large(&x_chunks, e0, jk, &mut qp);

        // compress qp into fw
        let mut fw = 0.0;
        for &qp_i in qp[..=jz].iter().rev() {
            fw = fw * f64::exp2i_fast(-24) + (qp_i as f64);
        }

        let y = scalbn(fw, qe.into()).set_sign(ih != 0);

        if x.is_sign_negative() {
            (n.wrapping_neg() & 3, -y)
        } else {
            (n & 3, y)
        }
    }
}

fn sin_core(x: f64) -> f32 {
    // GENERATE: sin_poly f64 4
    const K3: f64 = f64::from_bits(0xBFC5555554CAD521); // -1.6666666641473518e-1
    const K5: f64 = f64::from_bits(0x3F81111088D013ED); // 8.333329367824016e-3
    const K7: f64 = f64::from_bits(0xBF2A00F9613422C1); // -1.9839328943073965e-4
    const K9: f64 = f64::from_bits(0x3EC6CD684FD95971); // 2.7182546762746612e-6

    let x2 = x * x;
    let x3 = x2 * x;

    let r = x + horner!(x3, x2, [K3, K5, K7, K9]);

    r as f32
}

fn cos_core(x: f64) -> f32 {
    // GENERATE: cos_poly f64 4
    const K4: f64 = f64::from_bits(0x3FA5555554F45728); // 4.166666662255941e-2
    const K6: f64 = f64::from_bits(0xBF56C16B88F1DDCE); // -1.3888883729843318e-3
    const K8: f64 = f64::from_bits(0x3EFA01119AC045D9); // 2.4799513695999496e-5
    const K10: f64 = f64::from_bits(0xBE924297AC04A5CD); // -2.720970933173516e-7

    let x2 = x * x;
    let x4 = x2 * x2;

    let r = 1.0 - 0.5 * x2 + horner!(x4, x2, [K4, K6, K8, K10]);

    r as f32
}

fn tan_core(x: f64, inv: bool) -> f32 {
    // GENERATE: tan_poly f64 5
    const K3: f64 = f64::from_bits(0x3FD555555FAA85C5); // 3.333333429564756e-1
    const K5: f64 = f64::from_bits(0x3FC11109BE04155C); // 1.3333246019480327e-1
    const K7: f64 = f64::from_bits(0x3FABA516B9DAB93A); // 5.399390238504727e-2
    const K9: f64 = f64::from_bits(0x3F960F62619306C0); // 2.1543061452354406e-2
    const K11: f64 = f64::from_bits(0x3F85F0BCE827722E); // 1.071307738760038e-2

    // let y = 0.5 * x
    // tan(x) = 2 * tan(y) / (1 - tan(y)^2)

    let y = 0.5 * x;
    let y2 = y * y;
    let y3 = y2 * y;

    // t1 = tan(y)
    let t1 = y + horner!(y3, y2, [K3, K5, K7, K9, K11]);

    // t2 = 2 * tan(y)
    let t2 = 2.0 * t1;

    // t3 = 1 - tan(y)^2
    let t3 = 1.0 - t1 * t1;

    let r = if inv {
        // -1/tan(x) = -t3 / t2
        -(t3 / t2)
    } else {
        // tan(x) = t2 / t3
        t2 / t3
    };

    r as f32
}
