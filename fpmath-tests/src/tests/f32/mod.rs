mod cbrt;
mod exp;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod round;
mod sqrt;
mod trigonometric_deg;
mod trigonometric_pi;
mod trigonometric_rad;

#[derive(Debug)]
pub(crate) struct RefResult {
    exp: i32,
    hi: f32,
    lo: f32,
}

impl RefResult {
    #[cfg(test)]
    pub(crate) fn from_f64(value: f64) -> Self {
        use rustc_apfloat::ieee::{Double, Single};
        use rustc_apfloat::Float as _;
        use rustc_apfloat::FloatConvert as _;

        // Use soft-float to avoid X87 compiler bugs
        let mut tmp = Double::from_bits(value.to_bits().into());
        if tmp.is_nan() {
            Self {
                exp: 0,
                hi: f32::NAN,
                lo: f32::NAN,
            }
        } else if tmp > Single::largest().convert(&mut false).value {
            Self {
                exp: 0,
                hi: f32::INFINITY,
                lo: f32::INFINITY,
            }
        } else if tmp < -Single::largest().convert(&mut false).value {
            Self {
                exp: 0,
                hi: f32::NEG_INFINITY,
                lo: f32::NEG_INFINITY,
            }
        } else if tmp.is_zero() {
            Self {
                exp: 0,
                hi: 0.0,
                lo: 0.0,
            }
        } else {
            let e_real = tmp.ilogb() + 1;
            let e_sat = e_real.max(-126);
            tmp = tmp.scalbn(-e_sat);

            let hi: Single = tmp
                .convert_r(rustc_apfloat::Round::TowardZero, &mut false)
                .value;
            tmp = (tmp - hi.convert(&mut false).value).value;
            let lo: Single = tmp.convert(&mut false).value;

            Self {
                exp: e_sat,
                hi: f32::from_bits(hi.to_bits() as _),
                lo: f32::from_bits(lo.to_bits() as _),
            }
        }
    }

    #[cfg(test)]
    pub(crate) fn calc_error(&self, actual: f32) -> f32 {
        // Use MIN/MAX instead of infinity because with x87 there can be
        // non-infinity values greater than MAX/less than MIN.
        if actual.is_nan() || self.hi.is_nan() {
            f32::NAN
        } else if actual > f32::MAX {
            if self.hi > f32::MAX {
                0.0
            } else {
                f32::INFINITY
            }
        } else if actual < f32::MIN {
            if self.hi < f32::MIN {
                0.0
            } else {
                f32::INFINITY
            }
        } else if self.hi.abs() > f32::MAX {
            f32::INFINITY
        } else {
            let scaled = fpmath::scalbn(actual, -self.exp);
            let abs_err = ((scaled - self.hi) - self.lo).abs();
            fpmath::scalbn(abs_err, 24)
        }
    }
}

pub(crate) fn mkfloat(m: u32, e: i16, s: bool) -> f32 {
    let m = m >> (32 - 23);
    let e = u32::from((e + 127) as u16) << 23;
    let s = u32::from(s) << 31;
    f32::from_bits(m | e | s)
}

fn select_threshold(actual: f32, normal_th: f32, subnormal_th: f32) -> f32 {
    if actual == 0.0 || actual.is_subnormal() {
        subnormal_th
    } else {
        normal_th
    }
}

// Workaround X87 compiler bugs
fn purify(x: f32) -> f32 {
    std::hint::black_box(x)
}
