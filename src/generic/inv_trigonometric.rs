use crate::traits::{CastInto as _, Float, Int as _};

pub(crate) trait InvTrigonometric: Float {
    const PI: Self;
    const FRAC_PI_2: Self;

    fn asin_finite(x: Self) -> Self;

    fn acos_finite(x: Self) -> Self;

    fn atan_finite(x: Self) -> Self;

    fn atan2_finite(y: Self, x: Self) -> Self;

    fn asind_finite(x: Self) -> Self;

    fn acosd_finite(x: Self) -> Self;

    fn atand_finite(x: Self) -> Self;

    fn atan2d_finite(y: Self, x: Self) -> Self;

    fn asinpi_finite(x: Self) -> Self;

    fn acospi_finite(x: Self) -> Self;

    fn atanpi_finite(x: Self) -> Self;

    fn atan2pi_finite(y: Self, x: Self) -> Self;
}

pub(crate) fn asin<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asin(±1) = ±π/2
        F::FRAC_PI_2.copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero, asin(x) ~= x
        // also handles asin(-0) = -0
        x
    } else {
        F::asin_finite(x)
    }
}

pub(crate) fn acos<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.is_sign_negative() {
            // acos(-1) = π
            F::PI
        } else {
            // acos(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acos(x) ~= π/2
        F::FRAC_PI_2
    } else {
        F::acos_finite(x)
    }
}

pub(crate) fn atan<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atan(±inf) = ±π/2
            F::FRAC_PI_2.copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else if e <= F::MANT_BITS.into() {
        // subnormal or tiny, atan(x) ~= x
        // also handles atan(-0) = -0
        x
    } else {
        F::atan_finite(x)
    }
}

pub(crate) fn atan2<F: InvTrigonometric>(y: F, x: F) -> F {
    if x.is_nan() || y.is_nan() {
        // atan2(NaN, x) = NaN
        // atan2(y, NaN) = NaN
        return F::NAN;
    }

    let (ny, nx) = if y.raw_exp() <= F::MANT_BITS.into() || x.raw_exp() <= F::MANT_BITS.into() {
        // Convert possible subnormals to normals. This can overflow to infinity,
        // which `atan2_de_inf` handles below.
        let scale = F::exp2i_fast((F::MANT_BITS * 2 + 1).cast_into());
        (y * scale, x * scale)
    } else {
        (y, x)
    };
    let (ny, nx) = atan2_de_inf(ny, nx);

    let nxexp = nx.raw_exp();
    let nyexp = ny.raw_exp();
    if nyexp == F::RawExp::ZERO {
        if nx.is_sign_negative() {
            // atan2(±0, -0) = ±π
            // atan2(±0, x < 0) = ±π
            F::PI.copysign(ny)
        } else {
            // atan2(±0, +0) = ±0
            // atan2(±0, x > 0) = ±0
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // atan2(y < 0, ±0) = -π/2
        // atan2(y > 0, ±0) = +π/2
        F::FRAC_PI_2.copysign(ny)
    } else {
        F::atan2_finite(ny, nx)
    }
}

pub(crate) fn asind<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asind(±1) = ±90
        F::cast_from(90u8).copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else {
        F::asind_finite(x)
    }
}

pub(crate) fn acosd<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.is_sign_negative() {
            // acosd(-1) = 180
            F::cast_from(180u8)
        } else {
            // acosd(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acosd(x) ~= 90
        F::cast_from(90u8)
    } else {
        F::acosd_finite(x)
    }
}

pub(crate) fn atand<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atand(±inf) = ±90
            F::cast_from(90u8).copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else {
        F::atand_finite(x)
    }
}

pub(crate) fn atan2d<F: InvTrigonometric>(y: F, x: F) -> F {
    if x.is_nan() || y.is_nan() {
        // atan2d(NaN, x) = NaN
        // atan2d(y, NaN) = NaN
        return F::NAN;
    }

    let (ny, nx) = if y.raw_exp() <= F::MANT_BITS.into() || x.raw_exp() <= F::MANT_BITS.into() {
        // Convert possible subnormals to normals. This can overflow to infinity,
        // which `atan2_de_inf` handles below.
        let scale = F::exp2i_fast((F::MANT_BITS * 2 + 1).cast_into());
        (y * scale, x * scale)
    } else {
        (y, x)
    };
    let (ny, nx) = atan2_de_inf(ny, nx);

    let nxexp = nx.raw_exp();
    let nyexp = ny.raw_exp();
    if nyexp == F::RawExp::ZERO {
        if nx.is_sign_negative() {
            // atan2d(±0, -0) = ±180
            // atan2d(±0, x < 0) = ±180
            F::cast_from(180u8).copysign(ny)
        } else {
            // atan2d(±0, +0) = ±0
            // atan2d(±0, x > 0) = ±0
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // atan2d(y < 0, ±0) = -90
        // atan2d(y > 0, ±0) = +90
        F::cast_from(90u8).copysign(ny)
    } else {
        F::atan2d_finite(ny, nx)
    }
}

pub(crate) fn asinpi<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asinpi(±1) = ±0.5
        F::HALF.copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else {
        F::asinpi_finite(x)
    }
}

pub(crate) fn acospi<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.is_sign_negative() {
            // acospi(-1) = 1
            F::ONE
        } else {
            // acospi(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acospi(x) ~= 0.5
        F::HALF
    } else {
        F::acospi_finite(x)
    }
}

pub(crate) fn atanpi<F: InvTrigonometric>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atanpi(±inf) = ±0.5
            F::HALF.copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else {
        F::atanpi_finite(x)
    }
}

pub(crate) fn atan2pi<F: InvTrigonometric>(y: F, x: F) -> F {
    if x.is_nan() || y.is_nan() {
        // atan2pi(NaN, x) = NaN
        // atan2pi(y, NaN) = NaN
        return F::NAN;
    }

    let (ny, nx) = if y.raw_exp() <= F::MANT_BITS.into() || x.raw_exp() <= F::MANT_BITS.into() {
        // Convert possible subnormals to normals. This can overflow to infinity,
        // which `atan2_de_inf` handles below.
        let scale = F::exp2i_fast((F::MANT_BITS * 2 + 1).cast_into());
        (y * scale, x * scale)
    } else {
        (y, x)
    };
    let (ny, nx) = atan2_de_inf(ny, nx);

    let nxexp = nx.raw_exp();
    let nyexp = ny.raw_exp();
    if nyexp == F::RawExp::ZERO {
        if nx.is_sign_negative() {
            // atan2pi(±0, -0) = ±1
            // atan2pi(±0, x < 0) = ±1
            F::ONE.copysign(ny)
        } else {
            // atan2pi(±0, +0) = ±0
            // atan2pi(±0, x > 0) = ±0
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // atan2pi(y < 0, ±0) = -0.5
        // atan2pi(y > 0, ±0) = +0.5
        F::HALF.copysign(ny)
    } else {
        F::atan2pi_finite(ny, nx)
    }
}

#[inline]
fn atan2_de_inf<F: Float>(y: F, x: F) -> (F, F) {
    let y_inf = y.is_infinite();
    let x_inf = x.is_infinite();
    if y_inf || x_inf {
        let y = (if y_inf { F::ONE } else { F::ZERO }).copysign(y);
        let x = (if x_inf { F::ONE } else { F::ZERO }).copysign(x);
        (y, x)
    } else {
        (y, x)
    }
}

#[cfg(test)]
mod tests {
    use super::InvTrigonometric;
    use crate::FloatMath;

    fn test_asin<F: InvTrigonometric + FloatMath>() {
        use crate::asin;

        let f = F::parse;

        assert_is_nan!(asin(F::NAN));
        assert_is_nan!(asin(f("1.5")));
        assert_is_nan!(asin(f("-1.5")));
        assert_is_nan!(asin(F::INFINITY));
        assert_is_nan!(asin(F::NEG_INFINITY));
        assert_total_eq!(asin(F::ZERO), F::ZERO);
        assert_total_eq!(asin(-F::ZERO), -F::ZERO);
        assert_total_eq!(asin(F::ONE), F::FRAC_PI_2);
        assert_total_eq!(asin(-F::ONE), -F::FRAC_PI_2);
    }

    fn test_acos<F: InvTrigonometric + FloatMath>() {
        use crate::acos;

        let f = F::parse;

        assert_is_nan!(acos(F::NAN));
        assert_is_nan!(acos(f("1.5")));
        assert_is_nan!(acos(f("-1.5")));
        assert_is_nan!(acos(F::INFINITY));
        assert_is_nan!(acos(F::NEG_INFINITY));
        assert_total_eq!(acos(F::ZERO), F::FRAC_PI_2);
        assert_total_eq!(acos(-F::ZERO), F::FRAC_PI_2);
        assert_total_eq!(acos(F::ONE), F::ZERO);
        assert_total_eq!(acos(-F::ONE), F::PI);
    }

    fn test_atan<F: InvTrigonometric + FloatMath>() {
        use crate::atan;

        assert_is_nan!(atan(F::NAN));
        assert_total_eq!(atan(F::INFINITY), F::FRAC_PI_2);
        assert_total_eq!(atan(F::NEG_INFINITY), -F::FRAC_PI_2);
        assert_total_eq!(atan(F::ZERO), F::ZERO);
        assert_total_eq!(atan(-F::ZERO), -F::ZERO);
    }

    fn test_atan2<F: InvTrigonometric + FloatMath>() {
        use crate::{atan2, scalbn};

        assert_is_nan!(atan2(F::NAN, F::ONE));
        assert_is_nan!(atan2(F::NAN, F::ZERO));
        assert_is_nan!(atan2(F::NAN, F::INFINITY));
        assert_is_nan!(atan2(F::NAN, F::NAN));
        assert_is_nan!(atan2(F::INFINITY, F::NAN));
        assert_is_nan!(atan2(F::ZERO, F::NAN));
        assert_is_nan!(atan2(F::ONE, F::NAN));
        assert_total_eq!(atan2(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::ONE), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::ONE), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, -F::ZERO), F::PI);
        assert_total_eq!(atan2(-F::ZERO, -F::ZERO), -F::PI);
        assert_total_eq!(atan2(F::ZERO, -F::ONE), F::PI);
        assert_total_eq!(atan2(-F::ZERO, -F::ONE), -F::PI);
        assert_total_eq!(atan2(F::INFINITY, F::ZERO), F::FRAC_PI_2);
        assert_total_eq!(atan2(F::INFINITY, -F::ZERO), F::FRAC_PI_2);
        assert_total_eq!(atan2(F::INFINITY, F::ONE), F::FRAC_PI_2);
        assert_total_eq!(atan2(F::INFINITY, -F::ONE), F::FRAC_PI_2);
        assert_total_eq!(atan2(F::NEG_INFINITY, F::ZERO), -F::FRAC_PI_2);
        assert_total_eq!(atan2(F::NEG_INFINITY, -F::ZERO), -F::FRAC_PI_2);
        assert_total_eq!(atan2(F::NEG_INFINITY, F::ONE), -F::FRAC_PI_2);
        assert_total_eq!(atan2(F::NEG_INFINITY, -F::ONE), -F::FRAC_PI_2);
        assert_total_eq!(atan2(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::ONE, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::ONE, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::NEG_INFINITY), F::PI);
        assert_total_eq!(atan2(-F::ZERO, F::NEG_INFINITY), -F::PI);
        assert_total_eq!(atan2(F::ONE, F::NEG_INFINITY), F::PI);
        assert_total_eq!(atan2(-F::ONE, F::NEG_INFINITY), -F::PI);

        let small = scalbn(F::ONE, F::MIN_NORMAL_EXP.into() / 2);
        assert_total_eq!(atan2(small, F::LARGEST), F::ZERO);
        assert_total_eq!(atan2(-small, F::LARGEST), -F::ZERO);
        assert_total_eq!(atan2(small, -F::LARGEST), F::PI);
        assert_total_eq!(atan2(-small, -F::LARGEST), -F::PI);
    }

    fn test_asind<F: InvTrigonometric + FloatMath>() {
        use crate::asind;

        let f = F::parse;

        assert_is_nan!(asind(F::NAN));
        assert_is_nan!(asind(f("1.5")));
        assert_is_nan!(asind(f("-1.5")));
        assert_is_nan!(asind(F::INFINITY));
        assert_is_nan!(asind(F::NEG_INFINITY));
        assert_total_eq!(asind(F::ZERO), F::ZERO);
        assert_total_eq!(asind(-F::ZERO), -F::ZERO);
        assert_total_eq!(asind(F::ONE), f("90"));
        assert_total_eq!(asind(-F::ONE), f("-90"));
    }

    fn test_acosd<F: InvTrigonometric + FloatMath>() {
        use crate::acosd;

        let f = F::parse;

        assert_is_nan!(acosd(F::NAN));
        assert_is_nan!(acosd(f("1.5")));
        assert_is_nan!(acosd(f("-1.5")));
        assert_is_nan!(acosd(F::INFINITY));
        assert_is_nan!(acosd(F::NEG_INFINITY));
        assert_total_eq!(acosd(F::ZERO), f("90"));
        assert_total_eq!(acosd(-F::ZERO), f("90"));
        assert_total_eq!(acosd(F::ONE), F::ZERO);
        assert_total_eq!(acosd(-F::ONE), f("180"));
    }

    fn test_atand<F: InvTrigonometric + FloatMath>() {
        use crate::atand;

        let f = F::parse;

        assert_is_nan!(atand(F::NAN));
        assert_total_eq!(atand(F::INFINITY), f("90"));
        assert_total_eq!(atand(F::NEG_INFINITY), f("-90"));
        assert_total_eq!(atand(F::ZERO), F::ZERO);
        assert_total_eq!(atand(-F::ZERO), -F::ZERO);
    }

    fn test_atan2d<F: InvTrigonometric + FloatMath>() {
        use crate::{atan2d, scalbn};

        let f = F::parse;

        assert_is_nan!(atan2d(F::NAN, F::ONE));
        assert_is_nan!(atan2d(F::NAN, F::ZERO));
        assert_is_nan!(atan2d(F::NAN, F::INFINITY));
        assert_is_nan!(atan2d(F::NAN, F::NAN));
        assert_is_nan!(atan2d(F::INFINITY, F::NAN));
        assert_is_nan!(atan2d(F::ZERO, F::NAN));
        assert_is_nan!(atan2d(F::ONE, F::NAN));
        assert_total_eq!(atan2d(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::ONE), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::ONE), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, -F::ZERO), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::ZERO), f("-180"));
        assert_total_eq!(atan2d(F::ZERO, -F::ONE), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::ONE), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, F::ONE), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::ONE), f("90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, -F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::ONE), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, -F::ONE), f("-90"));
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ONE, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ONE, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::NEG_INFINITY), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, F::NEG_INFINITY), f("-180"));
        assert_total_eq!(atan2d(F::ONE, F::NEG_INFINITY), f("180"));
        assert_total_eq!(atan2d(-F::ONE, F::NEG_INFINITY), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::INFINITY), f("45"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::INFINITY), f("-45"));
        assert_total_eq!(atan2d(F::INFINITY, F::NEG_INFINITY), f("135"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::NEG_INFINITY), f("-135"));

        let small = scalbn(F::ONE, F::MIN_NORMAL_EXP.into() / 2);
        assert_total_eq!(atan2d(small, F::LARGEST), F::ZERO);
        assert_total_eq!(atan2d(-small, F::LARGEST), -F::ZERO);
        assert_total_eq!(atan2d(small, -F::LARGEST), f("180"));
        assert_total_eq!(atan2d(-small, -F::LARGEST), f("-180"));
    }

    fn test_asinpi<F: InvTrigonometric + FloatMath>() {
        use crate::asinpi;

        let f = F::parse;

        assert_is_nan!(asinpi(F::NAN));
        assert_is_nan!(asinpi(f("1.5")));
        assert_is_nan!(asinpi(f("-1.5")));
        assert_is_nan!(asinpi(F::INFINITY));
        assert_is_nan!(asinpi(F::NEG_INFINITY));
        assert_total_eq!(asinpi(F::ZERO), F::ZERO);
        assert_total_eq!(asinpi(-F::ZERO), -F::ZERO);
        assert_total_eq!(asinpi(F::ONE), F::HALF);
        assert_total_eq!(asinpi(-F::ONE), -F::HALF);
    }

    fn test_acospi<F: InvTrigonometric + FloatMath>() {
        use crate::acospi;

        let f = F::parse;

        assert_is_nan!(acospi(F::NAN));
        assert_is_nan!(acospi(f("1.5")));
        assert_is_nan!(acospi(f("-1.5")));
        assert_is_nan!(acospi(F::INFINITY));
        assert_is_nan!(acospi(F::NEG_INFINITY));
        assert_total_eq!(acospi(F::ZERO), F::HALF);
        assert_total_eq!(acospi(-F::ZERO), F::HALF);
        assert_total_eq!(acospi(F::ONE), F::ZERO);
        assert_total_eq!(acospi(-F::ONE), F::ONE);
    }

    fn test_atanpi<F: InvTrigonometric + FloatMath>() {
        use crate::atanpi;

        assert_is_nan!(atanpi(F::NAN));
        assert_total_eq!(atanpi(F::INFINITY), F::HALF);
        assert_total_eq!(atanpi(F::NEG_INFINITY), -F::HALF);
        assert_total_eq!(atanpi(F::ZERO), F::ZERO);
        assert_total_eq!(atanpi(-F::ZERO), -F::ZERO);
    }

    fn test_atan2pi<F: InvTrigonometric + FloatMath>() {
        use crate::{atan2pi, scalbn};

        let f = F::parse;

        assert_is_nan!(atan2pi(F::NAN, F::ONE));
        assert_is_nan!(atan2pi(F::NAN, F::ZERO));
        assert_is_nan!(atan2pi(F::NAN, F::INFINITY));
        assert_is_nan!(atan2pi(F::NAN, F::NAN));
        assert_is_nan!(atan2pi(F::INFINITY, F::NAN));
        assert_is_nan!(atan2pi(F::ZERO, F::NAN));
        assert_is_nan!(atan2pi(F::ONE, F::NAN));
        assert_total_eq!(atan2pi(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::ONE), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::ONE), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, -F::ZERO), F::ONE);
        assert_total_eq!(atan2pi(-F::ZERO, -F::ZERO), -F::ONE);
        assert_total_eq!(atan2pi(F::ZERO, -F::ONE), F::ONE);
        assert_total_eq!(atan2pi(-F::ZERO, -F::ONE), -F::ONE);
        assert_total_eq!(atan2pi(F::INFINITY, F::ZERO), F::HALF);
        assert_total_eq!(atan2pi(F::INFINITY, -F::ZERO), F::HALF);
        assert_total_eq!(atan2pi(F::INFINITY, F::ONE), F::HALF);
        assert_total_eq!(atan2pi(F::INFINITY, -F::ONE), F::HALF);
        assert_total_eq!(atan2pi(F::NEG_INFINITY, F::ZERO), -F::HALF);
        assert_total_eq!(atan2pi(F::NEG_INFINITY, -F::ZERO), -F::HALF);
        assert_total_eq!(atan2pi(F::NEG_INFINITY, F::ONE), -F::HALF);
        assert_total_eq!(atan2pi(F::NEG_INFINITY, -F::ONE), -F::HALF);
        assert_total_eq!(atan2pi(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::ONE, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::ONE, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::NEG_INFINITY), F::ONE);
        assert_total_eq!(atan2pi(-F::ZERO, F::NEG_INFINITY), -F::ONE);
        assert_total_eq!(atan2pi(F::ONE, F::NEG_INFINITY), F::ONE);
        assert_total_eq!(atan2pi(-F::ONE, F::NEG_INFINITY), -F::ONE);
        assert_total_eq!(atan2pi(F::INFINITY, F::INFINITY), f("0.25"));
        assert_total_eq!(atan2pi(F::NEG_INFINITY, F::INFINITY), f("-0.25"));
        assert_total_eq!(atan2pi(F::INFINITY, F::NEG_INFINITY), f("0.75"));
        assert_total_eq!(atan2pi(F::NEG_INFINITY, F::NEG_INFINITY), f("-0.75"));

        let small = scalbn(F::ONE, F::MIN_NORMAL_EXP.into() / 2);
        assert_total_eq!(atan2pi(small, F::LARGEST), F::ZERO);
        assert_total_eq!(atan2pi(-small, F::LARGEST), -F::ZERO);
        assert_total_eq!(atan2pi(small, -F::LARGEST), F::ONE);
        assert_total_eq!(atan2pi(-small, -F::LARGEST), -F::ONE);
    }

    #[test]
    fn test_f32() {
        test_asin::<f32>();
        test_acos::<f32>();
        test_atan::<f32>();
        test_atan2::<f32>();
        test_asind::<f32>();
        test_acosd::<f32>();
        test_atand::<f32>();
        test_atan2d::<f32>();
        test_asinpi::<f32>();
        test_acospi::<f32>();
        test_atanpi::<f32>();
        test_atan2pi::<f32>();
    }

    #[test]
    fn test_f64() {
        test_asin::<f64>();
        test_acos::<f64>();
        test_atan::<f64>();
        test_atan2::<f64>();
        test_asind::<f64>();
        test_acosd::<f64>();
        test_atand::<f64>();
        test_atan2d::<f64>();
        test_asinpi::<f64>();
        test_acospi::<f64>();
        test_atanpi::<f64>();
        test_atan2pi::<f64>();
    }
}
