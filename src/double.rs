use crate::traits::Float;

/// A double-float.
///
/// `hi` and `lo` might overlap partially.
#[derive(Copy, Clone, Debug)]
pub(crate) struct Double<F: Float> {
    hi: F,
    lo: F,
}

impl<F: Float> From<F> for Double<F> {
    #[inline]
    fn from(value: F) -> Self {
        Self {
            hi: value,
            lo: F::ZERO,
        }
    }
}

impl<F: Float> Double<F> {
    pub(crate) const ZERO: Self = Self {
        hi: F::ZERO,
        lo: F::ZERO,
    };

    #[inline]
    pub(crate) const fn new(hi: F, lo: F) -> Self {
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn hi(self) -> F {
        self.hi
    }

    #[inline]
    pub(crate) fn lo(self) -> F {
        self.lo
    }

    #[inline]
    pub(crate) fn to_single(self) -> F {
        self.hi + self.lo
    }

    #[inline]
    pub(crate) fn to_semi(self) -> SemiDouble<F> {
        let hi = self.hi.split_hi();
        let lo = (self.hi - hi) + self.lo;
        SemiDouble { hi, lo }
    }

    #[inline]
    pub(crate) fn normalize(self) -> Self {
        let hi = (self.hi + self.lo).purify();
        let lo = (self.hi - hi) + self.lo;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn abs(self) -> Self {
        if self.hi.sign() { -self } else { self }
    }

    #[inline]
    pub(crate) fn qadd1(self, rhs: F) -> Self {
        let hi = (self.hi + rhs).purify();
        let lo = ((self.hi - hi) + rhs) + self.lo;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn qradd1(self, rhs: F) -> Self {
        let hi = (self.hi + rhs).purify();
        let lo = ((rhs - hi) + self.hi) + self.lo;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn qadd2(self, rhs: Self) -> Self {
        let hi = (self.hi + rhs.hi).purify();
        let lo = ((self.hi - hi) + rhs.hi) + (self.lo + rhs.lo);
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn ladd(self, rhs_lo: F) -> Self {
        Self {
            hi: self.hi,
            lo: self.lo + rhs_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd11(lhs: F, rhs: F) -> Self {
        let hi = (lhs + rhs).purify();
        let lo = (lhs - hi) + rhs;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn new_add11(lhs: F, rhs: F) -> Self {
        let hi = (lhs + rhs).purify();
        let t1 = (hi - lhs).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs - t1).purify();
        let lo = (lhs - t2).purify() + t3;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn qsub1(self, rhs: F) -> Self {
        let hi = (self.hi - rhs).purify();
        let lo = ((self.hi - hi) - rhs) + self.lo;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn qsub2(self, rhs: Self) -> Self {
        let hi = (self.hi - rhs.hi).purify();
        let lo = ((self.hi - hi) - rhs.hi) + (self.lo - rhs.lo);
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn lsub(self, rhs_lo: F) -> Self {
        Self {
            hi: self.hi,
            lo: self.lo - rhs_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qsub11(lhs: F, rhs: F) -> Self {
        let hi = (lhs - rhs).purify();
        let lo = (lhs - hi) - rhs;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn new_sub11(lhs: F, rhs: F) -> Self {
        let hi = (lhs - rhs).purify();
        let t1 = (hi - lhs).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs + t1).purify();
        let lo = (lhs - t2).purify() - t3;
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn pmul1(self, rhs: F) -> Self {
        Self {
            hi: self.hi * rhs,
            lo: self.lo * rhs,
        }
    }

    #[inline]
    pub(crate) fn new_div11(lhs: F, rhs: F) -> Self {
        let (lhs_hi, lhs_lo) = lhs.split_hi_lo();
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let rhs_inv = F::ONE / rhs;
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = (lhs * rhs_inv).purify();
        let res_lo = -res_hi
            + lhs_hi * rhs_inv_hi
            + lhs_hi * rhs_inv_lo
            + lhs_lo * rhs_inv
            + res_hi * (F::ONE - rhs_hi * rhs_inv_hi - rhs_hi * rhs_inv_lo - rhs_lo * rhs_inv);

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn square(self) -> Self {
        self * self
    }

    #[inline]
    pub(crate) fn new_recip(rhs: F) -> Self {
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let rhs_inv = F::ONE / rhs;
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = rhs_inv.purify();
        let res_lo = -res_hi
            + rhs_inv_hi
            + rhs_inv_lo
            + res_hi * (F::ONE - rhs_hi * rhs_inv_hi - rhs_hi * rhs_inv_lo - rhs_lo * rhs_inv);

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn recip(self) -> Self {
        let rhs = self;
        let (rhs_hihi, rhs_hilo) = rhs.hi.split_hi_lo();
        let rhs_inv = (F::ONE / rhs.hi).purify();
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = rhs_inv;
        let res_lo = -res_hi
            + rhs_inv_hi
            + rhs_inv_lo
            + res_hi
                * (F::ONE
                    - rhs_hihi * rhs_inv_hi
                    - rhs_hihi * rhs_inv_lo
                    - rhs_hilo * rhs_inv_hi
                    - rhs_hilo * rhs_inv_lo);
        let res_lo = res_lo - res_hi * rhs.lo * rhs_inv;

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Neg for Double<F> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl<F: Float> core::ops::Add<F> for Double<F> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: F) -> Self {
        let hi = (self.hi + rhs).purify();
        let t1 = (hi - self.hi).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs - t1).purify();
        let lo = ((self.hi - t2).purify() + t3) + self.lo;
        Self { hi, lo }
    }
}

impl<F: Float> core::ops::Add for Double<F> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        let hi = (self.hi + rhs.hi).purify();
        let t1 = (hi - self.hi).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs.hi - t1).purify();
        let lo = ((self.hi - t2).purify() + t3) + (self.lo + rhs.lo);
        Self { hi, lo }
    }
}

impl<F: Float> core::ops::Sub<F> for Double<F> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: F) -> Self {
        let hi = (self.hi - rhs).purify();
        let t1 = (hi - self.hi).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs + t1).purify();
        let lo = ((self.hi - t2).purify() - t3) + self.lo;
        Self { hi, lo }
    }
}

impl<F: Float> core::ops::Sub for Double<F> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let hi = (self.hi - rhs.hi).purify();
        let t1 = (hi - self.hi).purify();
        let t2 = (hi - t1).purify();
        let t3 = (rhs.hi + t1).purify();
        let lo = ((self.hi - t2).purify() - t3) + (self.lo - rhs.lo);
        Self { hi, lo }
    }
}

impl<F: Float> core::ops::Mul for Double<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let lhs = self;
        let (lhs_hihi, lhs_hilo) = lhs.hi.split_hi_lo();
        let (rhs_hihi, rhs_hilo) = rhs.hi.split_hi_lo();

        let res_hi = (lhs.hi * rhs.hi).purify();
        let res_lo = lhs_hihi * rhs_hihi - res_hi
            + lhs_hilo * rhs_hihi
            + lhs_hihi * rhs_hilo
            + lhs_hilo * rhs_hilo
            + lhs.hi * rhs.lo
            + lhs.lo * rhs.hi;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Mul<F> for Double<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: F) -> Self {
        let lhs = self;
        let (lhs_hihi, lhs_hilo) = lhs.hi.split_hi_lo();
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let res_hi = (lhs.hi * rhs).purify();
        let res_lo = lhs_hihi * rhs_hi - res_hi
            + lhs_hilo * rhs_hi
            + lhs_hihi * rhs_lo
            + lhs_hilo * rhs_lo
            + lhs.lo * rhs;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Div for Double<F> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        let lhs = self;
        let (lhs_hihi, lhs_hilo) = lhs.hi.split_hi_lo();
        let (rhs_hihi, rhs_hilo) = rhs.hi.split_hi_lo();
        let rhs_inv = F::ONE / rhs.hi;
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = (lhs.hi * rhs_inv).purify();
        let res_lo = -res_hi
            + lhs_hihi * rhs_inv_hi
            + lhs_hihi * rhs_inv_lo
            + lhs_hilo * rhs_inv_hi
            + lhs_hilo * rhs_inv_lo
            + res_hi
                * (F::ONE
                    - rhs_hihi * rhs_inv_hi
                    - rhs_hihi * rhs_inv_lo
                    - rhs_hilo * rhs_inv_hi
                    - rhs_hilo * rhs_inv_lo);
        let res_lo = res_lo + (lhs.lo - res_hi * rhs.lo) * rhs_inv;

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

/// A semi-double-float.
///
/// The lower half of bits of `hi` are zero.
#[derive(Copy, Clone, Debug)]
pub(crate) struct SemiDouble<F: Float> {
    hi: F,
    lo: F,
}

impl<F: Float> SemiDouble<F> {
    pub(crate) const ONE: Self = Self {
        hi: F::ONE,
        lo: F::ZERO,
    };

    #[inline]
    pub(crate) fn new(value: F) -> Self {
        let (hi, lo) = value.split_hi_lo();
        Self { hi, lo }
    }

    #[inline]
    pub(crate) const fn with_parts(hi: F, lo: F) -> Self {
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn hi(self) -> F {
        self.hi
    }

    #[inline]
    pub(crate) fn to_single(self) -> F {
        self.hi + self.lo
    }

    #[inline]
    pub(crate) fn to_double(self) -> Double<F> {
        Double {
            hi: self.hi,
            lo: self.lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd11(lhs: F, rhs: F) -> Self {
        let res_hi = (lhs + rhs).split_hi();
        let res_lo = (lhs - res_hi) + rhs;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd21(lhs: Double<F>, rhs: F) -> Self {
        let res_hi = (lhs.hi + rhs).split_hi();
        let res_lo = ((lhs.hi - res_hi) + rhs) + lhs.lo;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd12(lhs: F, rhs: Double<F>) -> Self {
        let res_hi = (lhs + rhs.hi).split_hi();
        let res_lo = ((lhs - res_hi) + rhs.hi) + rhs.lo;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd22(lhs: Double<F>, rhs: Double<F>) -> Self {
        let res_hi = (lhs.hi + rhs.hi).split_hi();
        let res_lo = ((lhs.hi - res_hi) + rhs.hi) + (lhs.lo + rhs.lo);

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qsub11(lhs: F, rhs: F) -> Self {
        let res_hi = (lhs - rhs).split_hi();
        let res_lo = (lhs - res_hi) - rhs;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qsub12(lhs: F, rhs: Double<F>) -> Self {
        let res_hi = (lhs - rhs.hi).split_hi();
        let res_lo = ((lhs - res_hi) - rhs.hi) - rhs.lo;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn pmul1(self, rhs: F) -> Self {
        Self {
            hi: self.hi * rhs,
            lo: self.lo * rhs,
        }
    }

    #[inline]
    pub(crate) fn square(self) -> Double<F> {
        let res_hi = self.hi * self.hi;
        let res_lo = F::TWO * self.hi * self.lo + self.lo * self.lo;

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Neg for SemiDouble<F> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl<F: Float> core::ops::Mul for SemiDouble<F> {
    type Output = Double<F>;

    #[inline]
    fn mul(self, rhs: Self) -> Double<F> {
        let lhs = self;

        let res_hi = lhs.hi * rhs.hi;
        let res_lo = lhs.hi * rhs.lo + lhs.lo * rhs.hi + lhs.lo * rhs.lo;

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Mul<F> for SemiDouble<F> {
    type Output = Double<F>;

    #[inline]
    fn mul(self, rhs: F) -> Double<F> {
        let lhs = self;
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let res_hi = lhs.hi * rhs_hi;
        let res_lo = lhs.hi * rhs_lo + lhs.lo * rhs;

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Div for SemiDouble<F> {
    type Output = Double<F>;

    #[inline]
    fn div(self, rhs: Self) -> Double<F> {
        let lhs = self;
        let rhs_inv = F::ONE / (rhs.hi + rhs.lo).purify();
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = ((lhs.hi + lhs.lo) * rhs_inv).purify();
        let res_lo = -res_hi
            + lhs.hi * rhs_inv_hi
            + lhs.hi * rhs_inv_lo
            + lhs.lo * rhs_inv
            + res_hi * (F::ONE - rhs.hi * rhs_inv_hi - rhs.hi * rhs_inv_lo - rhs.lo * rhs_inv);

        Double {
            hi: res_hi,
            lo: res_lo,
        }
    }
}
