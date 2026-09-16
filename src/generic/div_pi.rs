use crate::double::SemiDouble;
use crate::traits::FloatConsts;

pub(crate) trait DivPi: FloatConsts {
    const FRAC_1_PI_EX: SemiDouble<Self>;
}
