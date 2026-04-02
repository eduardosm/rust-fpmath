use crate::double::SemiDouble;
use crate::traits::FloatConsts;

pub(crate) trait DivPi: FloatConsts {
    fn frac_1_pi_ex() -> SemiDouble<Self>;
}
