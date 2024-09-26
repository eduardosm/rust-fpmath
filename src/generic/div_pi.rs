use crate::double::SemiDouble;
use crate::traits::{FloatConsts, Like};

pub(crate) trait DivPi<L = Like<Self>>: FloatConsts {
    fn frac_1_pi_ex() -> SemiDouble<Self>;
}
