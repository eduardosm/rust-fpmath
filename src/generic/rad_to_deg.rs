use crate::double::SemiDouble;
use crate::traits::Float;

pub(crate) trait RadToDeg: Float {
    fn rad_to_deg() -> Self;
    fn rad_to_deg_ex() -> SemiDouble<Self>;
}
