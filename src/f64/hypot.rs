use super::f64x2::F64x2;
use crate::traits::Float as _;

impl crate::generic::Hypot for f64 {
    #[inline]
    fn hypot_finite(x: Self, y: Self) -> Self {
        let absx = x.abs();
        let absy = y.abs();
        let (min, max) = if absx < absy {
            (absx, absy)
        } else {
            (absy, absx)
        };
        if max == 0.0 {
            return 0.0;
        }

        let max_exp = max.exponent();
        let edelta = if max_exp >= 500 || max_exp <= -500 {
            (-max_exp).clamp(-1022, 1022)
        } else {
            0
        };

        let min = min * f64::exp2i_fast(edelta);
        let max = max * f64::exp2i_fast(edelta);

        (F64x2::square1(min) + F64x2::square1(max))
            .sqrt()
            .scalbn_to_f64((-edelta).into())
    }
}
