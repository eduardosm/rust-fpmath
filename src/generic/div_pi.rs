use crate::traits::{Float, Like};

pub(crate) trait DivPi<L = Like<Self>>: Float {
    fn frac_1_pi() -> Self;
    fn frac_1_pi_hi() -> Self;
    fn frac_1_pi_lo() -> Self;
}
