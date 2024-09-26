use crate::traits::Float;

/// A denormalized double-float.
///
/// `hi` and `lo` might overlap partially.
#[derive(Copy, Clone, Debug)]
pub(crate) struct DenormDouble<F: Float> {
    hi: F,
    lo: F,
}

impl<F: Float> DenormDouble<F> {
    #[inline]
    pub(crate) fn new(hi: F, lo: F) -> Self {
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
    pub(crate) fn to_norm(self) -> NormDouble<F> {
        let hi = (self.hi + self.lo).purify();
        let lo = (self.hi - hi) + self.lo;
        NormDouble { hi, lo }
    }

    #[inline]
    pub(crate) fn to_semi(self) -> SemiDouble<F> {
        let hi = self.hi.split_hi();
        let lo = (self.hi - hi) + self.lo;
        SemiDouble { hi, lo }
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
        let lo = if lhs > rhs {
            (lhs - hi) + rhs
        } else {
            (rhs - hi) + lhs
        };
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

        let rhs_inv = F::one() / rhs;
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = (lhs * rhs_inv).purify();
        let res_lo = -res_hi
            + lhs_hi * rhs_inv_hi
            + lhs_hi * rhs_inv_lo
            + lhs_lo * rhs_inv
            + res_hi * (F::one() - rhs_hi * rhs_inv_hi - rhs_hi * rhs_inv_lo - rhs_lo * rhs_inv);

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_recip(rhs: F) -> Self {
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let rhs_inv = F::one() / rhs;
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = rhs_inv.purify();
        let res_lo = -res_hi
            + rhs_inv_hi
            + rhs_inv_lo
            + res_hi * (F::one() - rhs_hi * rhs_inv_hi - rhs_hi * rhs_inv_lo - rhs_lo * rhs_inv);

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Neg for DenormDouble<F> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl<F: Float> core::ops::Add<F> for DenormDouble<F> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: F) -> Self {
        let hi = (self.hi + rhs).purify();
        let lo = if self.hi > rhs {
            ((self.hi - hi) + rhs) + self.lo
        } else {
            ((rhs - hi) + self.hi) + self.lo
        };
        Self { hi, lo }
    }
}

/// A normalized double-float.
///
/// `hi` and `lo` should not overlap.
#[derive(Copy, Clone, Debug)]
pub(crate) struct NormDouble<F: Float> {
    hi: F,
    lo: F,
}

impl<F: Float> NormDouble<F> {
    #[inline]
    pub(crate) fn with_parts(hi: F, lo: F) -> Self {
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
    pub(crate) fn to_denorm(self) -> DenormDouble<F> {
        DenormDouble {
            hi: self.hi,
            lo: self.lo,
        }
    }

    #[inline]
    pub(crate) fn qadd2(self, rhs: Self) -> DenormDouble<F> {
        let hi = (self.hi + rhs.hi).purify();
        let lo = ((self.hi - hi) + rhs.hi) + (self.lo + rhs.lo);
        DenormDouble { hi, lo }
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
    #[inline]
    pub(crate) fn new(value: F) -> Self {
        let (hi, lo) = value.split_hi_lo();
        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn one() -> Self {
        Self {
            hi: F::one(),
            lo: F::ZERO,
        }
    }

    #[inline]
    pub(crate) fn with_parts(hi: F, lo: F) -> Self {
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
    pub(crate) fn to_denorm(self) -> DenormDouble<F> {
        DenormDouble {
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
    pub(crate) fn new_qadd21(lhs: DenormDouble<F>, rhs: F) -> Self {
        let res_hi = (lhs.hi + rhs).split_hi();
        let res_lo = ((lhs.hi - res_hi) + rhs) + lhs.lo;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd12(lhs: F, rhs: DenormDouble<F>) -> Self {
        let res_hi = (lhs + rhs.hi).split_hi();
        let res_lo = ((lhs - res_hi) + rhs.hi) + rhs.lo;

        Self {
            hi: res_hi,
            lo: res_lo,
        }
    }

    #[inline]
    pub(crate) fn new_qadd22(lhs: DenormDouble<F>, rhs: DenormDouble<F>) -> Self {
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
    pub(crate) fn new_qsub12(lhs: F, rhs: DenormDouble<F>) -> Self {
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
    pub(crate) fn square(self) -> DenormDouble<F> {
        let res_hi = self.hi * self.hi;
        let res_lo = F::two() * self.hi * self.lo + self.lo * self.lo;

        DenormDouble {
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
    type Output = DenormDouble<F>;

    #[inline]
    fn mul(self, rhs: Self) -> DenormDouble<F> {
        let lhs = self;

        let res_hi = lhs.hi * rhs.hi;
        let res_lo = lhs.hi * rhs.lo + lhs.lo * rhs.hi + lhs.lo * rhs.lo;

        DenormDouble {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Mul<F> for SemiDouble<F> {
    type Output = DenormDouble<F>;

    #[inline]
    fn mul(self, rhs: F) -> DenormDouble<F> {
        let lhs = self;
        let (rhs_hi, rhs_lo) = rhs.split_hi_lo();

        let res_hi = lhs.hi * rhs_hi;
        let res_lo = lhs.hi * rhs_lo + lhs.lo * rhs;

        DenormDouble {
            hi: res_hi,
            lo: res_lo,
        }
    }
}

impl<F: Float> core::ops::Div for SemiDouble<F> {
    type Output = DenormDouble<F>;

    #[inline]
    fn div(self, rhs: Self) -> DenormDouble<F> {
        let lhs = self;
        let rhs_inv = F::one() / (rhs.hi + rhs.lo).purify();
        let (rhs_inv_hi, rhs_inv_lo) = rhs_inv.split_hi_lo();

        let res_hi = ((lhs.hi + lhs.lo) * rhs_inv).purify();
        let res_lo = -res_hi
            + lhs.hi * rhs_inv_hi
            + lhs.hi * rhs_inv_lo
            + lhs.lo * rhs_inv
            + res_hi * (F::one() - rhs.hi * rhs_inv_hi - rhs.hi * rhs_inv_lo - rhs.lo * rhs_inv);

        DenormDouble {
            hi: res_hi,
            lo: res_lo,
        }
    }
}
