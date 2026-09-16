use crate::double::SemiDouble;
use crate::traits::Float;

pub(crate) trait RadToDeg: Float {
    const RAD_TO_DEG: Self;
    const RAD_TO_DEG_EX: SemiDouble<Self>;
}
