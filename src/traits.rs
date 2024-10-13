pub(crate) trait CastFrom<T> {
    fn cast_from(value: T) -> Self;
}

pub(crate) trait CastInto<T> {
    fn cast_into(self) -> T;
}

impl<T, U: CastFrom<T>> CastInto<U> for T {
    #[inline]
    fn cast_into(self) -> U {
        U::cast_from(self)
    }
}

pub(crate) trait Int:
    'static
    + Copy
    + Ord
    + core::fmt::Debug
    + From<bool>
    + CastInto<i8>
    + CastInto<u8>
    + CastInto<i16>
    + CastInto<u16>
    + CastInto<i32>
    + CastInto<u32>
    + CastInto<i64>
    + CastInto<u64>
    + CastFrom<i8>
    + CastFrom<u8>
    + CastFrom<i16>
    + CastFrom<u16>
    + CastFrom<i32>
    + CastFrom<u32>
    + CastFrom<i64>
    + CastFrom<u64>
    + TryFrom<u8>
    + TryInto<i32>
    + core::ops::Add<Self, Output = Self>
    + core::ops::Sub<Self, Output = Self>
    + core::ops::Mul<Self, Output = Self>
    + core::ops::Div<Self, Output = Self>
    + core::ops::Rem<Self, Output = Self>
    + core::ops::AddAssign<Self>
    + core::ops::SubAssign<Self>
    + core::ops::Not<Output = Self>
    + core::ops::BitAnd<Self, Output = Self>
    + core::ops::BitOr<Self, Output = Self>
    + core::ops::BitXor<Self, Output = Self>
    + core::ops::Shl<Self, Output = Self>
    + core::ops::Shr<Self, Output = Self>
    + core::ops::Shl<u8, Output = Self>
    + core::ops::Shr<u8, Output = Self>
    + core::ops::ShlAssign<u8>
    + core::ops::ShrAssign<u8>
{
    const ZERO: Self;
    const ONE: Self;
    const TWO: Self;

    const MAX: Self;
}

#[allow(dead_code)] // https://github.com/rust-lang/rust/issues/128839
pub(crate) trait UInt: Int + From<u8> {}

#[allow(dead_code)] // https://github.com/rust-lang/rust/issues/128839
pub(crate) trait SInt: Int + From<i8> + core::ops::Neg<Output = Self> {}

pub(crate) trait Float:
    'static
    + Copy
    + PartialOrd
    + CastFrom<u8>
    + CastFrom<i16>
    + CastFrom<i32>
    + CastFrom<u32>
    + CastFrom<i64>
    + CastFrom<u64>
    + core::fmt::Debug
    + core::fmt::Display
    + core::ops::Neg<Output = Self>
    + core::ops::Add<Self, Output = Self>
    + core::ops::Sub<Self, Output = Self>
    + core::ops::Mul<Self, Output = Self>
    + core::ops::Div<Self, Output = Self>
{
    // Hack to avoid "conflicting implementations of trait"
    type Like;

    type Raw: UInt
        + From<Self::RawExp>
        + From<u16>
        + CastInto<Self>
        + core::ops::Shl<Self::RawExp, Output = Self::Raw>
        + core::ops::Shl<Self::Exp, Output = Self::Raw>
        + core::ops::Shr<Self::RawExp, Output = Self::Raw>
        + core::ops::Shr<Self::Exp, Output = Self::Raw>;

    type RawExp: UInt + CastFrom<Self::Raw>;

    type Exp: SInt + CastInto<Self> + Into<i32>;

    const BITS: u8;
    const MANT_BITS: u8;
    const EXP_BITS: u8;

    const SIGN_MASK: Self::Raw;
    const EXP_MASK: Self::Raw;
    const MANT_MASK: Self::Raw;

    const EXP_OFFSET: Self::RawExp;
    const MAX_RAW_EXP: Self::RawExp;

    const MIN_NORMAL_EXP: Self::Exp;
    const MAX_EXP: Self::Exp;

    const INFINITY: Self;
    fn neg_infinity() -> Self;
    const NAN: Self;

    const ZERO: Self;
    fn half() -> Self;
    fn one() -> Self;
    fn two() -> Self;

    #[cfg(test)]
    fn largest() -> Self;

    /// Workarounds rustc/LLVM bugs
    fn purify(self) -> Self;

    fn to_raw(self) -> Self::Raw;

    fn from_raw(raw: Self::Raw) -> Self;

    fn raw_exp_to_exp(e: Self::RawExp) -> Self::Exp;

    fn exp_to_raw_exp(e: Self::Exp) -> Self::RawExp;

    #[inline]
    fn sign(self) -> bool {
        (self.to_raw() & Self::SIGN_MASK) != Self::Raw::ZERO
    }

    #[inline]
    fn raw_exp(self) -> Self::RawExp {
        ((self.to_raw() & Self::EXP_MASK) >> Self::MANT_BITS).cast_into()
    }

    // Not named `exp` to avoid confussion with the exp function
    #[inline]
    fn exponent(self) -> Self::Exp {
        Self::raw_exp_to_exp(self.raw_exp())
    }

    #[inline]
    fn normalize_arg(self) -> (Self, Self::Exp) {
        if self.raw_exp() == Self::RawExp::ZERO {
            // convert possible subnormal to normal
            let escale = Self::Exp::cast_from(Self::MANT_BITS);
            (self * Self::exp2i_fast(escale), -escale)
        } else {
            (self, Self::Exp::ZERO)
        }
    }

    #[inline]
    fn raw_mant(self) -> Self::Raw {
        self.to_raw() & Self::MANT_MASK
    }

    #[inline]
    fn mant(self) -> Self::Raw {
        (self.to_raw() & Self::MANT_MASK) | (Self::MANT_MASK + Self::Raw::ONE)
    }

    #[cfg(test)]
    fn is_nan(self) -> bool;

    #[inline]
    fn abs(self) -> Self {
        Self::from_raw(self.to_raw() & !Self::SIGN_MASK)
    }

    #[inline]
    fn copysign(self, y: Self) -> Self {
        Self::from_raw((self.to_raw() & !Self::SIGN_MASK) | (y.to_raw() & Self::SIGN_MASK))
    }

    #[inline]
    fn set_sign(self, s: bool) -> Self {
        Self::from_raw(
            (self.to_raw() & !Self::SIGN_MASK)
                | (Self::Raw::from(s) << (Self::MANT_BITS + Self::EXP_BITS)),
        )
    }

    #[inline]
    fn set_raw_exp(self, e: Self::RawExp) -> Self {
        Self::from_raw((self.to_raw() & !Self::EXP_MASK) | (Self::Raw::from(e) << Self::MANT_BITS))
    }

    /// Returns a float with the mantissa and sign of `x`
    /// and the exponent `e`
    ///
    /// `MIN_NORMAL_EXP <= e <= MAX_EXP`
    #[inline]
    fn set_exp(self, e: Self::Exp) -> Self {
        self.set_raw_exp(Self::exp_to_raw_exp(e))
    }

    #[inline]
    fn exp2i_fast(x: Self::Exp) -> Self {
        Self::one().set_exp(x)
    }

    #[inline]
    fn split_hi(self) -> Self {
        Self::from_raw(self.to_raw() & (Self::Raw::MAX << ((Self::MANT_BITS + 2) / 2)))
    }

    #[inline]
    fn split_hi_lo(self) -> (Self, Self) {
        let x = self.purify();
        let hi = x.split_hi();
        let lo = x - hi;
        (hi, lo)
    }

    #[inline]
    fn norm_hi_lo_full(hi: Self, lo: Self) -> (Self, Self) {
        let lo = lo.purify();
        let hi2 = (hi + lo).purify();
        let lo2 = (hi - hi2) + lo;
        (hi2, lo2)
    }

    #[inline]
    fn norm_hi_lo_splitted(hi: Self, lo: Self) -> (Self, Self) {
        let lo = lo.purify();
        let hi2 = hi.split_hi();
        let lo2 = (hi - hi2) + lo;
        (hi2, lo2)
    }

    #[cfg(test)]
    fn parse(s: &str) -> Self;
}

// Hack to avoid "conflicting implementations of trait"
pub(crate) type Like<F> = <F as Float>::Like;

pub(crate) trait FloatConsts<L = Like<Self>>: Float {
    fn pi() -> Self;
    fn frac_pi_2() -> Self;
    fn frac_pi_4() -> Self;
    fn frac_2_pi() -> Self;
}
