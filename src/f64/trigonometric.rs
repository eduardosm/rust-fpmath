use super::f64x2::F64x2;
use crate::generic::{reduce_90_deg, reduce_half_revs, reduce_pi_2_large, round_fi};
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_PI_180 PI
const FRAC_PI_180: F64x2 = F64x2::from_bits(0x3F91DF46A2529D39, 0x3C15C1D8BECDD291); // 1.745329251994329576923690768489e-2
const PI: F64x2 = F64x2::from_bits(0x400921FB54442D18, 0x3CA1A62633145C07); // 3.141592653589793238462643383280e0

impl crate::generic::Trigonometric for f64 {
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
        if x.exponent() <= -40 {
            // Tiny, tan(x) ~= x
            x
        } else {
            let (n, y) = reduce_pi_2(x);
            tan_core(y, (n & 1) != 0)
        }
    }

    fn sind_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_90_deg(x);
            if y == 0.0 {
                match n {
                    0 => 0.0f64.copysign(x),
                    1 => 1.0,
                    2 => 0.0f64.copysign(x),
                    3 => -1.0,
                    _ => unreachable!(),
                }
            } else {
                let z = y * FRAC_PI_180;
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
            let z = y * FRAC_PI_180;
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
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            let sin = (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106);
            (sin, 1.0)
        } else {
            let (n, y) = reduce_90_deg(x);
            if y == 0.0 {
                match n {
                    0 => (0.0f64.copysign(x), 1.0),
                    1 => (1.0, 0.0),
                    2 => (0.0f64.copysign(x), -1.0),
                    3 => (-1.0, 0.0),
                    _ => unreachable!(),
                }
            } else {
                let z = y * FRAC_PI_180;
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
        if x.exponent() <= -40 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_90_deg(x);
            let inv = (n & 1) != 0;
            if y == 0.0 {
                if inv {
                    f64::INFINITY.set_sign(n == 3)
                } else {
                    0.0f64.set_sign(x.is_sign_negative() ^ (n == 2))
                }
            } else {
                let z = y * FRAC_PI_180;
                tan_core(z, inv)
            }
        }
    }

    fn sinpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_half_revs(x);
            if y == 0.0 {
                match n {
                    0 => 0.0f64.copysign(x),
                    1 => 1.0,
                    2 => 0.0f64.copysign(x),
                    3 => -1.0,
                    _ => unreachable!(),
                }
            } else {
                let z = y * PI;
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
            let z = y * PI;
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
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            let sin = (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106);
            (sin, 1.0)
        } else {
            let (n, y) = reduce_half_revs(x);
            if y == 0.0 {
                match n {
                    0 => (0.0f64.copysign(x), 1.0),
                    1 => (1.0, 0.0),
                    2 => (0.0f64.copysign(x), -1.0),
                    3 => (-1.0, 0.0),
                    _ => unreachable!(),
                }
            } else {
                let z = y * PI;
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
        if x.exponent() <= -40 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_half_revs(x);
            let inv = (n & 1) != 0;
            if y == 0.0 {
                if inv {
                    f64::INFINITY.set_sign(n == 3)
                } else {
                    0.0f64.set_sign(x.is_sign_negative() ^ (n == 2))
                }
            } else {
                let z = y * PI;
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
fn reduce_pi_2(x: f64) -> (u8, F64x2) {
    let x_orig = x;

    // GENERATE: consts f64 FRAC_PI_4 FRAC_2_PI
    const FRAC_PI_4: f64 = f64::from_bits(0x3FE921FB54442D18); // 7.853981633974483e-1
    const FRAC_2_PI: f64 = f64::from_bits(0x3FE45F306DC9C883); // 6.366197723675814e-1

    // GENERATE: reduce_pi_2::medium_consts f64 F64x2
    const FRAC_PI_2_P0: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
    const FRAC_PI_2_P0EX: F64x2 = F64x2::from_bits(0x3C91A62633145C07, 0xB91F1976B7ED8FBC); // 6.123233995736765886130329661375e-17

    if x.abs() <= FRAC_PI_4 {
        // reduction not needed
        return (0, F64x2::new1(x));
    }

    let p = x * FRAC_2_PI;
    if p.exponent() <= 40 {
        let (f_n, n) = round_fi(p);

        let x = F64x2::new1(x);
        let f_n = F64x2::new1(f_n);

        let r = x - f_n * FRAC_PI_2_P0;
        let y = r - f_n * FRAC_PI_2_P0EX;

        (n as u8 & 3, y)
    } else {
        let mant = x_orig.mant();
        let x_chunks = [
            ((mant >> 29) as u32),
            (((mant >> 5) & 0x00FF_FFFF) as u32),
            (((mant << 19) & 0x00FF_FFFF) as u32),
        ];
        let e0 = x_orig.exponent() - 23;
        let jk = 6;

        let mut qp: [u64; 20] = [0; 20];
        let (ih, jz, n, qe) = reduce_pi_2_large(&x_chunks, e0, jk, &mut qp);

        // compress qp into fw
        let mut fw = F64x2::ZERO;
        for &qp_i in qp[..=jz].iter().rev() {
            fw = fw.scalbn_fast(-24) + (qp_i as f64);
        }

        let mut y = fw.scalbn_medium(qe.into());
        if ih != 0 {
            y = -y;
        }

        if x.is_sign_negative() {
            (n.wrapping_neg() & 3, -y)
        } else {
            (n & 3, y)
        }
    }
}

fn sin_core(x: F64x2) -> f64 {
    // GENERATE: sin_poly F64x2 7
    const K3: F64x2 = F64x2::from_bits(0xBFC5555555555555, 0xBC64BDE0078E4CE5); // -1.666666666666666664100864650420e-1
    const K5: F64x2 = F64x2::from_bits(0x3F8111111111110B, 0xBBF70A79A932D4D6); // 8.333333333333322731278575972540e-3
    const K7: F64x2 = F64x2::from_bits(0xBF2A01A01A018ABF, 0xBBC46506AC48AFA9); // -1.984126984125502378318557385877e-4
    const K9: F64x2 = F64x2::from_bits(0x3EC71DE3A533670D, 0xBB6DC568DCE6F92A); // 2.755731921416719049540893847479e-6
    const K11: F64x2 = F64x2::from_bits(0xBE5AE6452944750E, 0xBAD98F8BC5562307); // -2.505210490541038385806980909121e-8
    const K13: F64x2 = F64x2::from_bits(0x3DE6120901C62D2B, 0x3A729751AC67BB96); // 1.605836582519096210535731054289e-10
    const K15: F64x2 = F64x2::from_bits(0xBD6AAA3F69241CA9, 0xB9E12DA6C6BE4CB0); // -7.578657366197662834739917002720e-13

    let x2 = x.square();
    let x3 = x2 * x;

    let r = x + horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15]);

    r.to_f64()
}

fn cos_core(x: F64x2) -> f64 {
    // GENERATE: cos_poly F64x2 7
    const K4: F64x2 = F64x2::from_bits(0x3FA5555555555555, 0x3C44D9203F8833C1); // 4.166666666666666661406280190441e-2
    const K6: F64x2 = F64x2::from_bits(0xBF56C16C16C16C10, 0x3BD3F7678DC32046); // -1.388888888888887407099068108923e-3
    const K8: F64x2 = F64x2::from_bits(0x3EFA01A01A018DCA, 0x3B60D7DD23F9F2B7); // 2.480158730157141811542179442969e-5
    const K10: F64x2 = F64x2::from_bits(0xBE927E4FB75FB51A, 0xBB38DE8625E89DFC); // -2.755731921534173911168928805425e-7
    const K12: F64x2 = F64x2::from_bits(0x3E21EED8CA138E70, 0xBAACFCD8F9337BF8); // 2.087675435832826491402066604875e-9
    const K14: F64x2 = F64x2::from_bits(0xBDA939335745AC35, 0xBA46A91D5C78C2C5); // -1.147029236111312457594409803024e-11
    const K16: F64x2 = F64x2::from_bits(0x3D2AAC465E13328E, 0xB9B88F48B3F536D2); // 4.738067489025742534441685245595e-14

    let x2 = x.square();
    let x4 = x2.square();

    let r = 1.0 - x2.halve() + horner!(x4, x2, [K4, K6, K8, K10, K12, K14, K16]);

    r.to_f64()
}

fn tan_core(x: F64x2, inv: bool) -> f64 {
    // GENERATE: tan_poly F64x2 10
    const K3: F64x2 = F64x2::from_bits(0x3FD5555555555555, 0xBC7241B68EA6A5B8); // 3.333333333333332989944602813442e-1
    const K5: F64x2 = F64x2::from_bits(0x3FC1111111111294, 0x3C6086AA5284E9A5); // 1.333333333333440800577507064898e-1
    const K7: F64x2 = F64x2::from_bits(0x3FABA1BA1B9F2ADD, 0x3C3676C8DE2CD542); // 5.396825396709002114689637485274e-2
    const K9: F64x2 = F64x2::from_bits(0x3F9664F4893E5463, 0xBC302B9BAB8E4BA3); // 2.186948859851546143654712216091e-2
    const K11: F64x2 = F64x2::from_bits(0x3F8226E3140A1BA4, 0xBBCE391F24D4BB4F); // 8.863233613068129111390625388824e-3
    const K13: F64x2 = F64x2::from_bits(0x3F6D6D50930D5CC6, 0xBC0395B3D378CE4A); // 3.592164394301630397048373454530e-3
    const K15: F64x2 = F64x2::from_bits(0x3F57D85EAE67AF08, 0x3BEA44AE98A3B5E6); // 1.455395185639698055293487535981e-3
    const K17: F64x2 = F64x2::from_bits(0x3F4371D50C944DF3, 0x3BE8FCACFF68BF96); // 5.934038428300744889519843402900e-4
    const K19: F64x2 = F64x2::from_bits(0x3F2D42FBDADE40E0, 0x3BB27926168B9117); // 2.232487144593403275133374376594e-4
    const K21: F64x2 = F64x2::from_bits(0x3F220D17A5340579, 0x3BBC1A97954EB410); // 1.377192844316634906109434973675e-4

    // let y = 0.5 * x
    // tan(x) = 2 * tan(y) / (1 - tan(y)^2)

    let y = x.halve();
    let y2 = y.square();
    let y3 = y2 * y;

    // t1 = tan(y)
    let t1 = y + horner!(y3, y2, [K3, K5, K7, K9, K11, K13, K15, K17, K19, K21]);

    // t2 = 2 * tan(y)
    let t2 = t1.twice();

    // t3 = 1 - tan(y)^2
    let t3 = 1.0 - t1.square();

    let r = if inv {
        // -1/tan(x) = -t3 / t2
        -(t3 / t2)
    } else {
        // tan(x) = t2 / t3
        t2 / t3
    };

    r.to_f64()
}
