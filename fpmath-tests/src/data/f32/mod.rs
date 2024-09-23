pub(crate) mod acosh;
pub(crate) mod asin_acos;
pub(crate) mod asinh;
pub(crate) mod atan;
pub(crate) mod atan2;
pub(crate) mod atanh;
pub(crate) mod cbrt;
pub(crate) mod exp;
pub(crate) mod hypot;
pub(crate) mod log;
pub(crate) mod log_1p;
pub(crate) mod pow;
pub(crate) mod powi;
pub(crate) mod sin_cos;
pub(crate) mod sind_cosd;
pub(crate) mod sinh_cosh;
pub(crate) mod sinpi_cospi;
pub(crate) mod sqrt;
pub(crate) mod tan;
pub(crate) mod tand;
pub(crate) mod tanh;
pub(crate) mod tanpi;

#[inline]
pub(crate) fn mkfloat(m: u32, e: i16, s: bool) -> f32 {
    let m = m >> (32 - 23);
    let e = u32::from((e + 127) as u16) << 23;
    let s = u32::from(s) << 31;
    f32::from_bits(m | e | s)
}

#[derive(Debug, bincode::Encode, bincode::Decode)]
pub(crate) struct RefResult {
    exp: i32,
    hi: f32,
    lo: f32,
}

impl RefResult {
    pub(crate) fn from_rug(mut tmp: rug::Float) -> Self {
        if tmp.is_nan() {
            Self {
                exp: 0,
                hi: f32::NAN,
                lo: f32::NAN,
            }
        } else if tmp > f32::MAX {
            Self {
                exp: 0,
                hi: f32::INFINITY,
                lo: f32::INFINITY,
            }
        } else if tmp < f32::MIN {
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
            let e_real = tmp.get_exp().unwrap();
            let e_sat = e_real.max(-126);
            tmp >>= e_sat;

            let hi = tmp.to_f32_round(rug::float::Round::Zero);
            tmp -= hi;
            let lo = tmp.to_f32();

            Self { exp: e_sat, hi, lo }
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

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct OneArgData {
    pub(crate) x: f32,
    pub(crate) expected: RefResult,
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct TwoArgData {
    pub(crate) x: f32,
    pub(crate) y: f32,
    pub(crate) expected: RefResult,
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct SinCosData {
    pub(crate) x: f32,
    pub(crate) expected_sin: RefResult,
    pub(crate) expected_cos: RefResult,
}
