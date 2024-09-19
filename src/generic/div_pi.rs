use crate::traits::{FloatConsts, Like};

pub(crate) trait DivPi<L = Like<Self>>: FloatConsts {
    fn frac_1_pi_hi() -> Self;
    fn frac_1_pi_lo() -> Self;
}
