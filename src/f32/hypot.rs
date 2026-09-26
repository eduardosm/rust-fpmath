impl crate::generic::Hypot for f32 {
    #[inline]
    fn hypot_finite(x: Self, y: Self) -> Self {
        let x = f64::from(x);
        let y = f64::from(y);
        crate::f64::fast_sqrt(x * x + y * y).0 as f32
    }
}
