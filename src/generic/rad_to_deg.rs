use crate::traits::{Float, Like};

pub(crate) trait RadToDeg<L = Like<Self>>: Float {
    fn rad_to_deg() -> Self;
    fn rad_to_deg_hi() -> Self;
    fn rad_to_deg_lo() -> Self;
}
