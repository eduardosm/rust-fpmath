#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]

use std::os::raw::{c_int, c_long, c_ulong};

#[allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    unreachable_pub,
    clippy::redundant_static_lifetimes
)]
mod ffi {
    include!(concat!(env!("OUT_DIR"), "/mpfr_ffi.rs"));
}

pub type Prec = ffi::mpfr_prec_t;

pub type RawMpfr = ffi::mpfr_t;

pub type Exp = ffi::mpfr_exp_t;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Rnd {
    N,
    Z,
    U,
    D,
    A,
    NA,
}

impl Rnd {
    #[inline]
    fn as_raw(self) -> ffi::mpfr_rnd_t {
        match self {
            Self::N => ffi::MPFR_RNDN,
            Self::Z => ffi::MPFR_RNDZ,
            Self::U => ffi::MPFR_RNDU,
            Self::D => ffi::MPFR_RNDD,
            Self::A => ffi::MPFR_RNDA,
            Self::NA => ffi::MPFR_RNDNA,
        }
    }
}

pub struct Mpfr {
    raw: ffi::mpfr_t,
}

impl Drop for Mpfr {
    #[inline]
    fn drop(&mut self) {
        unsafe {
            ffi::mpfr_clear(self.as_mut_ptr());
        }
    }
}

impl Default for Mpfr {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for Mpfr {
    #[inline]
    fn clone(&self) -> Self {
        let mut x = Self::new2(self.get_prec());
        x.set(self, Rnd::Z);
        x
    }
}

impl Mpfr {
    #[inline]
    pub fn new() -> Self {
        unsafe {
            let mut raw: ffi::mpfr_t = std::mem::zeroed();
            ffi::mpfr_init(raw.as_mut_ptr());
            Self::from_raw(raw)
        }
    }

    #[inline]
    pub fn new2(prec: Prec) -> Self {
        unsafe {
            let mut raw: ffi::mpfr_t = std::mem::zeroed();
            ffi::mpfr_init2(raw.as_mut_ptr(), prec);
            Self::from_raw(raw)
        }
    }

    #[inline]
    pub fn from_raw(raw: RawMpfr) -> Self {
        Self { raw }
    }

    #[inline]
    pub fn as_ptr(&self) -> ffi::mpfr_srcptr {
        self.raw.as_ptr()
    }

    #[inline]
    pub fn as_mut_ptr(&mut self) -> ffi::mpfr_ptr {
        self.raw.as_mut_ptr()
    }

    #[inline]
    pub fn set_prec(&mut self, prec: Prec) {
        unsafe {
            ffi::mpfr_set_prec(self.as_mut_ptr(), prec);
        }
    }

    #[inline]
    pub fn get_prec(&self) -> Prec {
        unsafe { ffi::mpfr_get_prec(self.as_ptr()) }
    }

    #[inline]
    pub fn is_nan(&self) -> bool {
        unsafe { ffi::mpfr_nan_p(self.as_ptr()) != 0 }
    }

    #[inline]
    pub fn is_inf(&self) -> bool {
        unsafe { ffi::mpfr_inf_p(self.as_ptr()) != 0 }
    }

    #[inline]
    pub fn is_regular(&self) -> bool {
        unsafe { ffi::mpfr_regular_p(self.as_ptr()) != 0 }
    }

    #[inline]
    pub fn set_zero(&mut self, sign: bool) {
        unsafe {
            ffi::mpfr_set_zero(self.as_mut_ptr(), if sign { -1 } else { 1 });
        }
    }

    #[inline]
    pub fn set_inf(&mut self, sign: bool) {
        unsafe {
            ffi::mpfr_set_inf(self.as_mut_ptr(), if sign { -1 } else { 1 });
        }
    }

    #[inline]
    pub fn const_pi(&mut self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_const_pi(self.as_mut_ptr(), rnd.as_raw());
        }
    }

    #[inline]
    pub fn set(&mut self, value: &Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_set(self.as_mut_ptr(), value.as_ptr(), rnd.as_raw());
        }
    }

    #[inline]
    pub fn set_si(&mut self, value: c_long, rnd: Rnd) {
        unsafe {
            ffi::mpfr_set_si(self.as_mut_ptr(), value, rnd.as_raw());
        }
    }

    #[inline]
    pub fn set_ui(&mut self, value: c_ulong, rnd: Rnd) {
        unsafe {
            ffi::mpfr_set_ui(self.as_mut_ptr(), value, rnd.as_raw());
        }
    }

    #[inline]
    pub fn set_f32(&mut self, value: f32, rnd: Rnd) {
        unsafe {
            ffi::mpfr_set_flt(self.as_mut_ptr(), value, rnd.as_raw());
        }
    }

    #[inline]
    pub fn set_f64(&mut self, value: f64, rnd: Rnd) {
        unsafe {
            ffi::mpfr_set_d(self.as_mut_ptr(), value, rnd.as_raw());
        }
    }

    #[inline]
    pub fn set_exp(&mut self, e: Exp) {
        unsafe {
            ffi::mpfr_set_exp(self.as_mut_ptr(), e);
        }
    }

    #[inline]
    pub fn get_ui(&self, rnd: Rnd) -> c_ulong {
        unsafe { ffi::mpfr_get_ui(self.as_ptr(), rnd.as_raw()) }
    }

    #[inline]
    pub fn get_f32(&self, rnd: Rnd) -> f32 {
        unsafe { ffi::mpfr_get_flt(self.as_ptr(), rnd.as_raw()) }
    }

    #[inline]
    pub fn get_f64(&self, rnd: Rnd) -> f64 {
        unsafe { ffi::mpfr_get_d(self.as_ptr(), rnd.as_raw()) }
    }

    #[inline]
    pub fn get_exp(&self) -> Exp {
        assert!(self.is_regular());
        unsafe { ffi::mpfr_get_exp(self.as_ptr()) }
    }

    #[inline]
    pub fn sgn(&self) -> c_int {
        unsafe { ffi::mpfr_sgn(self.as_ptr()) }
    }

    #[inline]
    pub fn rint(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_rint(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn round(&mut self, x: Option<&Self>) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_round(self.as_mut_ptr(), x);
        }
    }

    #[inline]
    pub fn trunc(&mut self, x: Option<&Self>) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_trunc(self.as_mut_ptr(), x);
        }
    }

    #[inline]
    pub fn ceil(&mut self, x: Option<&Self>) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_ceil(self.as_mut_ptr(), x);
        }
    }

    #[inline]
    pub fn floor(&mut self, x: Option<&Self>) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_floor(self.as_mut_ptr(), x);
        }
    }

    #[inline]
    pub fn cmp_f32(&self, y: f32) -> std::cmp::Ordering {
        unsafe { ffi::mpfr_cmp_d(self.as_ptr(), f64::from(y)).cmp(&0) }
    }

    #[inline]
    pub fn cmp_f64(&self, y: f64) -> std::cmp::Ordering {
        unsafe { ffi::mpfr_cmp_d(self.as_ptr(), y).cmp(&0) }
    }

    #[inline]
    pub fn sub_ui(&mut self, x: Option<&Self>, y: c_ulong, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_sub_ui(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn sub_f32(&mut self, x: Option<&Self>, y: f32, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_sub_d(self.as_mut_ptr(), x, f64::from(y), rnd.as_raw());
        }
    }

    #[inline]
    pub fn sub_f64(&mut self, x: Option<&Self>, y: f64, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_sub_d(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn mul(&mut self, x: Option<&Self>, y: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            let y = y.map_or_else(|| self.as_ptr(), |y| y.as_ptr());
            ffi::mpfr_mul(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn div(&mut self, x: Option<&Self>, y: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            let y = y.map_or_else(|| self.as_ptr(), |y| y.as_ptr());
            ffi::mpfr_div(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn fmod(&mut self, x: Option<&Self>, y: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            let y = y.map_or_else(|| self.as_ptr(), |y| y.as_ptr());
            ffi::mpfr_fmod(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn mul_ui(&mut self, x: Option<&Self>, y: c_ulong, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_mul_ui(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn mul_2ui(&mut self, x: Option<&Self>, y: c_ulong, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_mul_2ui(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn mul_f64(&mut self, x: Option<&Self>, y: f64, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_mul_d(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn div_f64(&mut self, x: Option<&Self>, y: f64, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_div_d(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn f64_div(&mut self, x: f64, y: Option<&Self>, rnd: Rnd) {
        unsafe {
            let y = y.map_or_else(|| self.as_ptr(), |y| y.as_ptr());
            ffi::mpfr_d_div(self.as_mut_ptr(), x, y, rnd.as_raw());
        }
    }

    #[inline]
    pub fn abs(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_abs(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn sqrt(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_sqrt(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn cbrt(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_cbrt(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn hypot(&mut self, x: &Self, y: &Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_hypot(self.as_mut_ptr(), x.as_ptr(), y.as_ptr(), rnd.as_raw());
        }
    }

    #[inline]
    pub fn log(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_log(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn log2(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_log2(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn log10(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_log10(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn log1p(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_log1p(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn exp(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_exp(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn exp2(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_exp2(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn exp10(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_exp10(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn expm1(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_expm1(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn pow(&mut self, x: &Self, y: &Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_pow(self.as_mut_ptr(), x.as_ptr(), y.as_ptr(), rnd.as_raw());
        }
    }

    #[inline]
    pub fn sin_cos(&self, sin_out: &mut Self, cos_out: &mut Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_sin_cos(
                sin_out.as_mut_ptr(),
                cos_out.as_mut_ptr(),
                self.as_ptr(),
                rnd.as_raw(),
            );
        }
    }

    #[inline]
    pub fn tan(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_tan(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn sinh_cosh(&self, sinh_out: &mut Self, cosh_out: &mut Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_sinh_cosh(
                sinh_out.as_mut_ptr(),
                cosh_out.as_mut_ptr(),
                self.as_ptr(),
                rnd.as_raw(),
            );
        }
    }

    #[inline]
    pub fn tanh(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_tanh(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn asin(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_asin(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn acos(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_acos(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn atan(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_atan(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn atan2(&mut self, y: &Self, x: &Self, rnd: Rnd) {
        unsafe {
            ffi::mpfr_atan2(self.as_mut_ptr(), y.as_ptr(), x.as_ptr(), rnd.as_raw());
        }
    }

    #[inline]
    pub fn asinh(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_asinh(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn acosh(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_acosh(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }

    #[inline]
    pub fn atanh(&mut self, x: Option<&Self>, rnd: Rnd) {
        unsafe {
            let x = x.map_or_else(|| self.as_ptr(), |x| x.as_ptr());
            ffi::mpfr_atanh(self.as_mut_ptr(), x, rnd.as_raw());
        }
    }
}
