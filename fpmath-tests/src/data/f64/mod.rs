pub(crate) mod acosh;
pub(crate) mod asin_acos;
pub(crate) mod asinh;
pub(crate) mod atan;
pub(crate) mod atan2;
pub(crate) mod atanh;
pub(crate) mod cbrt;
pub(crate) mod exp;
pub(crate) mod gamma;
pub(crate) mod hypot;
pub(crate) mod log;
pub(crate) mod log_1p;
pub(crate) mod pow;
pub(crate) mod powi;
pub(crate) mod sin_cos;
pub(crate) mod sind_cosd;
pub(crate) mod sinh_cosh;
pub(crate) mod sinpi_cospi;
pub(crate) mod tan;
pub(crate) mod tand;
pub(crate) mod tanh;
pub(crate) mod tanpi;

#[inline]
pub(crate) fn mkfloat(m: u64, e: i16, s: bool) -> f64 {
    let m = m >> (64 - 52);
    let e = u64::from((e + 1023) as u16) << 52;
    let s = u64::from(s) << 63;
    f64::from_bits(m | e | s)
}

const RUG_PREC: u32 = 53 + 20;

#[derive(Debug, bincode::Encode, bincode::Decode)]
pub(crate) struct RefResult {
    exp: i32,
    hi: f64,
    lo: f64,
}

impl RefResult {
    pub(crate) fn from_rug(mut tmp: rug::Float) -> Self {
        if tmp.is_nan() {
            Self {
                exp: 0,
                hi: f64::NAN,
                lo: f64::NAN,
            }
        } else if tmp > f64::MAX {
            Self {
                exp: 0,
                hi: f64::INFINITY,
                lo: f64::INFINITY,
            }
        } else if tmp < f64::MIN {
            Self {
                exp: 0,
                hi: f64::NEG_INFINITY,
                lo: f64::NEG_INFINITY,
            }
        } else if tmp.is_zero() {
            Self {
                exp: 0,
                hi: 0.0,
                lo: 0.0,
            }
        } else {
            let e_real = tmp.get_exp().unwrap();
            let e_sat = e_real.max(-1022);
            tmp >>= e_sat;

            let hi = tmp.to_f64_round(rug::float::Round::Zero);
            tmp -= hi;
            let lo = tmp.to_f64();

            Self { exp: e_sat, hi, lo }
        }
    }

    #[cfg(test)]
    pub(crate) fn calc_error(&self, actual: f64) -> f64 {
        // Use MIN/MAX instead of infinity because with x87 there can be
        // non-infinity values greater than MAX/less than MIN.
        if actual.is_nan() && self.hi.is_nan() {
            0.0
        } else if actual.is_nan() || self.hi.is_nan() {
            f64::NAN
        } else if actual > f64::MAX {
            if self.hi > f64::MAX {
                0.0
            } else {
                f64::INFINITY
            }
        } else if actual < f64::MIN {
            if self.hi < f64::MIN {
                0.0
            } else {
                f64::INFINITY
            }
        } else if self.hi.abs() > f64::MAX {
            f64::INFINITY
        } else {
            let scaled = fpmath::scalbn(actual, -self.exp);
            let abs_err = ((scaled - self.hi) - self.lo).abs();
            fpmath::scalbn(abs_err, 53)
        }
    }
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct OneArgData {
    pub(crate) x: f64,
    pub(crate) expected: RefResult,
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct TwoArgData {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) expected: RefResult,
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct SinCosData {
    pub(crate) x: f64,
    pub(crate) expected_sin: RefResult,
    pub(crate) expected_cos: RefResult,
}
