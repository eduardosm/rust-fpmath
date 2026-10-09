#![allow(clippy::suspicious_arithmetic_impl)]

use crate::traits::Float as _;

// Users of this type are responsible for avoiding infinities and NaNs
// (and operations that produce them).

#[derive(Copy, Clone, Debug)]
pub(super) struct F64x2 {
    hi: f64,
    lo: f64,
}

impl F64x2 {
    #[inline]
    pub(super) const fn new(hi: f64, lo: f64) -> Self {
        // This constructor is meant to be used in constants,
        // so the assert does not impact runtime performance.
        assert!(hi.is_finite() && lo.is_finite());
        let hi_exp = (hi.to_bits() >> 52) & 0x7FF;
        let half_ulp = if hi_exp >= 54 {
            f64::from_bits((hi_exp - 53) << 52) // normal
        } else if hi_exp >= 2 {
            f64::from_bits(1u64 << (hi_exp - 2)) // subnormal
        } else {
            0.0 // not representable, so `lo` must be zero
        };
        assert!(lo == 0.0 || lo.abs() <= half_ulp);
        Self { hi, lo }
    }

    #[inline]
    pub(super) const fn from_bits(hi: u64, lo: u64) -> Self {
        Self::new(f64::from_bits(hi), f64::from_bits(lo))
    }

    #[inline]
    pub(super) const fn new1(hi: f64) -> Self {
        Self { hi, lo: 0.0 }
    }

    #[inline]
    pub(super) fn hi(self) -> f64 {
        self.hi
    }

    #[inline]
    pub(super) fn lo(self) -> f64 {
        self.lo
    }

    #[inline]
    pub(super) fn to_f64(self) -> f64 {
        purify(self.hi + self.lo)
    }

    #[inline]
    pub(super) fn abs(self) -> Self {
        if self.hi.is_sign_negative() {
            Self {
                hi: -self.hi,
                lo: -self.lo,
            }
        } else {
            self
        }
    }

    #[inline]
    pub(super) fn add11(lhs: f64, rhs: f64) -> Self {
        let (hi, lo) = two_sum(lhs, rhs);
        Self { hi, lo }
    }

    /// Assumes `|lhs| >= |rhs|` or that `lhs + rhs` is exact.
    #[inline]
    pub(super) fn fast_add11(lhs: f64, rhs: f64) -> Self {
        let (hi, lo) = fast_two_sum(lhs, rhs);
        Self { hi, lo }
    }

    #[inline]
    pub(super) fn sub11(lhs: f64, rhs: f64) -> Self {
        let (hi, lo) = two_sub(lhs, rhs);
        Self { hi, lo }
    }

    #[inline]
    pub(super) fn mul11(lhs: f64, rhs: f64) -> Self {
        let (hi, lo) = two_prod(lhs, rhs);
        Self { hi, lo }
    }

    #[inline]
    pub(super) fn recip(self) -> Self {
        let rhs_inv = purify(self.hi.recip());
        let (p, pe) = two_prod(self.hi, rhs_inv);
        let rh = purify(purify(1.0 - p) - pe);
        let rl = -purify(self.lo * rhs_inv);
        let (eh, el) = two_sum(rh, rl);
        let e = Self { hi: eh, lo: el } * rhs_inv;
        e + rhs_inv
    }

    #[inline]
    pub(super) fn div11(lhs: f64, rhs: f64) -> Self {
        let q = purify(lhs / rhs);
        let (p, pe) = two_prod(q, rhs);
        let r = purify(purify(lhs - p) - pe);
        let (hi, lo) = fast_two_sum(q, purify(r / rhs));
        Self { hi, lo }
    }

    /// Like `scalbn(self.to_f64(), y)`, but rounding only once, so
    /// subnormals do not suffer double rounding.
    #[inline]
    pub(super) fn scalbn_to_f64(self, y: i32) -> f64 {
        let min_normal = crate::generic::scalbn(f64::MIN_POSITIVE, -y);
        if self.hi.abs() >= min_normal || min_normal.is_infinite() {
            // Result is normal or infinite
            crate::generic::scalbn(self.to_f64(), y)
        } else {
            // Result is subnormal, adjust for correct rounding.
            // Copysign is used to preserve the sign when the result
            // underflows to zero.
            let c = min_normal.copysign(self.hi);
            let r = purify((self + c).to_f64() - c).copysign(self.hi);
            crate::generic::scalbn(r, y)
        }
    }
}

/// A number `f = f1 + f2`, where `f1` has 26 significant bits, so its
/// products with other numbers of 26 bits are exact.
#[derive(Copy, Clone)]
pub(super) struct Split {
    f: f64,
    f1: f64,
    f2: f64,
}

impl Split {
    #[inline]
    pub(super) fn new(f: f64) -> Self {
        let f1 = f.split_hi();
        Self { f, f1, f2: f - f1 }
    }

    /// Returns `f1`.
    #[inline]
    pub(super) fn hi(self) -> f64 {
        self.f1
    }

    /// Returns `f2 = f - f1`.
    #[inline]
    pub(super) fn lo(self) -> f64 {
        self.f2
    }

    /// Returns `(hi, lo)` such that `hi + lo ~= f * (a_hi + a_lo)`, with a
    /// relative error of about 2^-76 + 2^-52 * |a_lo / a_hi|.
    ///
    /// `hi` has at most 52 bits and `|lo|` is at most about
    /// `(2^-24 + |a_lo / a_hi|) * |hi|`.
    #[inline]
    pub(super) fn mul(self, a_hi: f64, a_lo: f64) -> (f64, f64) {
        // a_hi = a1 + a2, where `a1` has 26 bits, so `a1 * f1` and `a1 * f2`
        // are exact. The rounding errors of `a2 * f` and the first sum (less
        // than 2^-24 * |f * a_hi|) add about 2^-76, and those of `a_lo * f`
        // and the second sum add about 2^-52 * |a_lo / a_hi|.
        let a1 = a_hi.split_hi();
        let a2 = a_hi - a1;
        (a1 * self.f1, (a1 * self.f2 + a2 * self.f) + a_lo * self.f)
    }

    /// Returns `(hi, lo)` such that `hi + lo ~= c + f * (a_hi + a_lo)`, where
    /// `|f * (a_hi + a_lo)| <= |c.hi()|`, with the product calculated by
    /// `mul`.
    ///
    /// When `c` and the product have the same sign, `|lo|` is at most about
    /// `(2^-24 + |a_lo / a_hi|) / 2 * |hi|`.
    #[inline]
    pub(super) fn mul_add(self, a_hi: f64, a_lo: f64, c: F64x2) -> (f64, f64) {
        let (p_hi, p_lo) = self.mul(a_hi, a_lo);
        let s = F64x2::fast_add11(c.hi(), p_hi);
        (s.hi(), s.lo() + (p_lo + c.lo()))
    }
}

impl core::ops::Neg for F64x2 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl core::ops::Add for F64x2 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        let (sh, sl) = two_sum(self.hi, rhs.hi);
        let (th, tl) = two_sum(self.lo, rhs.lo);
        let c = purify(sl + th);
        let (vh, vl) = fast_two_sum(sh, c);
        let w = purify(tl + vl);
        let (hi, lo) = fast_two_sum(vh, w);
        Self { hi, lo }
    }
}

impl core::ops::Add<f64> for F64x2 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: f64) -> Self {
        let (sh, sl) = two_sum(self.hi, rhs);
        let c = purify(sl + self.lo);
        let (hi, lo) = fast_two_sum(sh, c);
        Self { hi, lo }
    }
}

impl core::ops::Add<F64x2> for f64 {
    type Output = F64x2;

    #[inline]
    fn add(self, rhs: F64x2) -> F64x2 {
        rhs + self
    }
}

impl core::ops::AddAssign for F64x2 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl core::ops::AddAssign<f64> for F64x2 {
    #[inline]
    fn add_assign(&mut self, rhs: f64) {
        *self = *self + rhs;
    }
}

impl core::ops::Sub for F64x2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let (sh, sl) = two_sub(self.hi, rhs.hi);
        let (th, tl) = two_sub(self.lo, rhs.lo);
        let c = purify(sl + th);
        let (vh, vl) = fast_two_sum(sh, c);
        let w = purify(tl + vl);
        let (hi, lo) = fast_two_sum(vh, w);
        Self { hi, lo }
    }
}

impl core::ops::Sub<F64x2> for f64 {
    type Output = F64x2;

    #[inline]
    fn sub(self, rhs: F64x2) -> F64x2 {
        let (sh, sl) = two_sub(self, rhs.hi);
        let c = purify(sl - rhs.lo);
        let (hi, lo) = fast_two_sum(sh, c);
        F64x2 { hi, lo }
    }
}

impl core::ops::Sub<f64> for F64x2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: f64) -> Self {
        let (sh, sl) = two_sub(self.hi, rhs);
        let c = purify(sl + self.lo);
        let (hi, lo) = fast_two_sum(sh, c);
        Self { hi, lo }
    }
}

impl core::ops::SubAssign for F64x2 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl core::ops::SubAssign<f64> for F64x2 {
    #[inline]
    fn sub_assign(&mut self, rhs: f64) {
        *self = *self - rhs;
    }
}

impl core::ops::Mul for F64x2 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let (ch, cl) = two_prod(self.hi, rhs.hi);
        let thl = purify(self.hi * rhs.lo);
        let tlh = purify(self.lo * rhs.hi);
        let (hi, lo) = fast_two_sum(ch, purify(purify(thl + tlh) + cl));
        Self { hi, lo }
    }
}

impl core::ops::Mul<f64> for F64x2 {
    type Output = F64x2;

    #[inline]
    fn mul(self, rhs: f64) -> Self {
        let (ch, cl) = two_prod(self.hi, rhs);
        let (hi, lo) = fast_two_sum(ch, purify(purify(self.lo * rhs) + cl));
        Self { hi, lo }
    }
}

impl core::ops::Mul<F64x2> for f64 {
    type Output = F64x2;

    #[inline]
    fn mul(self, rhs: F64x2) -> F64x2 {
        rhs * self
    }
}

impl core::ops::MulAssign for F64x2 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl core::ops::MulAssign<f64> for F64x2 {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

impl core::ops::Div for F64x2 {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.recip()
    }
}

impl core::ops::DivAssign for F64x2 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = purify(a + b);
    let t = purify(s - a);
    let e = purify(purify(a - purify(s - t)) + purify(b - t));
    (s, e)
}

#[inline]
fn two_sub(a: f64, b: f64) -> (f64, f64) {
    let s = purify(a - b);
    let t = purify(s - a);
    let e = purify(purify(a - purify(s - t)) - purify(b + t));
    (s, e)
}

/// Assumes `|lhs| >= |rhs|` or that `lhs + rhs` is exact.
#[inline]
fn fast_two_sum(lhs: f64, rhs: f64) -> (f64, f64) {
    let s = purify(lhs + rhs);
    let e = purify(rhs - purify(s - lhs));
    (s, e)
}

#[inline]
fn two_prod(lhs: f64, rhs: f64) -> (f64, f64) {
    let (lhs_hi, lhs_lo) = crate::traits::Float::split_hi_lo(lhs);
    let (rhs_hi, rhs_lo) = crate::traits::Float::split_hi_lo(rhs);

    // `lhs_lo` and `rhs_lo` have up to 27 significant bits, so
    // `lhs_lo * rhs_lo` would need 54 bits to be exact. Split `rhs_lo`
    // again so the multiplication can be exact.
    let rhs_lo_hi = f64::from_bits(rhs_lo.to_bits() & (u64::MAX << 39));
    let rhs_lo_lo = purify(rhs_lo - rhs_lo_hi);

    let hh = purify(lhs_hi * rhs_hi);
    let hl = purify(lhs_hi * rhs_lo);
    let lh = purify(lhs_lo * rhs_hi);
    let llh = purify(lhs_lo * rhs_lo_hi);
    let lll = purify(lhs_lo * rhs_lo_lo);

    let rh = purify(lhs * rhs);
    let rl = purify(purify(purify(purify(purify(hh - rh) + hl) + lh) + llh) + lll);

    (rh, rl)
}

#[inline]
fn purify(x: f64) -> f64 {
    crate::traits::Float::purify(x)
}
