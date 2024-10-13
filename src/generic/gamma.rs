use super::exp::{exp_split, hi_lo_exp_inner_common};
use super::log::hi_lo_log_hi_lo_inner;
use super::scalbn;
use super::sin_cos::{hi_lo_cos_inner, hi_lo_sin_inner};
use super::{is_int, reduce_half_mul_pi, Exp, Log, ReduceHalfMulPi, SinCos};
use crate::double::{DenormDouble, NormDouble, SemiDouble};
use crate::traits::{Float, FloatConsts, Int as _, Like};

pub(crate) trait Gamma<L = Like<Self>>:
    FloatConsts + SinCos + ReduceHalfMulPi + Exp + Log
{
    fn lo_th() -> Self;
    fn hi_th() -> Self;

    fn th_1() -> Self;
    fn th_2() -> Self;
    fn th_3() -> Self;

    const POLY_OFF: u8;

    fn half_ln_2_pi() -> NormDouble<Self>;

    fn lgamma_poly_1(x: Self) -> (Self, Self, Self, Self);
    fn lgamma_poly_2(x: Self) -> (Self, Self, Self, Self);

    fn special_poly(x: Self) -> Self;
}

pub(crate) fn tgamma<F: Gamma>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // tgamma(±0) = ±inf
        F::INFINITY.copysign(x)
    } else if x >= F::hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if e == F::MAX_RAW_EXP {
        // tgamma(NaN or -inf) = NaN
        F::NAN
    } else if x.sign() && is_int(x) {
        // tgamma(neg integer) = NaN
        F::NAN
    } else if x < F::lo_th() {
        // -inf and negative integers are handled above
        F::ZERO
    } else {
        tgamma_inner(x)
    }
}

pub(crate) fn lgamma<F: Gamma>(x: F) -> (F, i8) {
    let e = x.raw_exp();
    let sign = x.sign();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // lgamma(0) = inf
        (F::INFINITY, if sign { -1 } else { 1 })
    } else if e == F::MAX_RAW_EXP {
        if !sign && x.raw_mant() == F::Raw::ZERO {
            // lgamma(inf) = inf
            (F::INFINITY, 1)
        } else {
            // lgamma(NaN or -inf) = NaN
            (F::NAN, 0)
        }
    } else if x.sign() && is_int(x) {
        // lgamma(neg integer) = inf
        (F::INFINITY, 0)
    } else if x == F::one() || x == F::two() {
        // lgamma(1 or 2) = 0
        // ensure positive zero
        (F::ZERO, 1)
    } else {
        lgamma_inner(x)
    }
}

fn tgamma_inner<F: Gamma>(x: F) -> F {
    let (y, s) = gamma_inner_common(x);

    let (k, r_hi, r_lo) = exp_split(y.hi());
    let r_lo = r_lo + y.lo();
    let exp_y = hi_lo_exp_inner_common(r_hi, r_lo);

    scalbn((s.to_semi() * exp_y.to_semi()).to_single(), k)
}

fn lgamma_inner<F: Gamma>(x: F) -> (F, i8) {
    let (y, s) = gamma_inner_common(x);

    if y.hi() == F::INFINITY {
        (F::INFINITY, 1)
    } else {
        let s = s.to_norm();
        let (sign, abs_s) = if s.hi().sign() { (-1, -s) } else { (1, s) };

        let log_s = hi_lo_log_hi_lo_inner(abs_s, F::Exp::ZERO);

        (log_s.qadd2(y).to_single(), sign)
    }
}

/// Returns `(y, s)` such as `Γ(x) = s * exp(y)`.
fn gamma_inner_common<F: Gamma>(x: F) -> (DenormDouble<F>, DenormDouble<F>) {
    // For x < 0.5, use gamma reflection formula:
    // Γ(x)*Γ(1-x) = π/sin(πx) => Γ(x) = π/(sin(πx)*Γ(1-x))
    let reflect = (x < F::half()).then(|| {
        let (n, z) = reduce_half_mul_pi(x);
        let sinpix = match n {
            0 => hi_lo_sin_inner(z),
            1 => hi_lo_cos_inner(z),
            2 => -hi_lo_sin_inner(z),
            3 => -hi_lo_cos_inner(z),
            _ => unreachable!(),
        };
        // π / sin(πx)
        F::pi_ex() / sinpix.to_semi()
    });
    // nx is always greater or equal to 0.5
    let nx = if reflect.is_some() {
        DenormDouble::new_sub11(F::one(), x)
    } else {
        DenormDouble::new(x, F::ZERO)
    };

    // Based on the algorithm used in SLEEF.

    if nx.hi() < F::th_2() {
        // For small values of `nx`, ln(Γ(nx)) is calculated using a polynomial.

        let nx = nx.hi();
        let y;
        let r;
        let k1;
        let k2;
        let k3;
        if nx < F::th_1() {
            y = nx - F::one();
            (r, k1, k2, k3) = F::lgamma_poly_1(y);
        } else {
            y = nx - F::two();
            (r, k1, k2, k3) = F::lgamma_poly_2(y);
        };
        // r = ln(Γ(nx))
        let r = finish_poly(y, r, k1, k2, k3);

        if let Some(reflect) = reflect {
            // -ln(Γ(1 - x)), π / sin(πx)
            (-r, reflect)
        } else {
            // ln(Γ(x)), 1
            (r, DenormDouble::one())
        }
    } else {
        // For larger values of `nx`:
        // t = nx or nx + POLY_OFF
        // Γ(nx) = (P(1 / t) / t + 1) * t^(t - 0.5) * e^(-t) * √(2π)
        // P is a polynomial.

        let low = nx.hi() < F::th_3();
        let t = if low {
            nx + F::cast_from(F::POLY_OFF)
        } else {
            nx
        };
        let tinv = F::one() / t.to_single();

        // p = P(1 / t) * (1 / t) + 1
        let p1 = F::special_poly(tinv);
        let p = SemiDouble::new(p1) * SemiDouble::new(tinv) + F::one();

        // r = (t - 0.5) * ln(t) - t + 0.5 * ln(2π)
        //   = t * (ln(t) - 1) - 0.5 * ln(t) + 0.5 * ln(2π)
        let log_t = hi_lo_log_hi_lo_inner(t.to_norm(), F::Exp::ZERO);
        let r = t.to_semi() * (log_t - F::one()).to_semi() - log_t.pmul1(F::half())
            + F::half_ln_2_pi().to_denorm();

        let s = if low {
            let mut den = nx;
            for i in 1..F::POLY_OFF {
                let nx_plus_i = nx + F::cast_from(i);
                den = (den * nx_plus_i).normalize();
            }
            if let Some(reflect) = reflect {
                (reflect * den) / p
            } else {
                p / den
            }
        } else if let Some(reflect) = reflect {
            reflect / p
        } else {
            p
        };

        (if reflect.is_some() { -r } else { r }, s)
    }
}

fn finish_poly<F: Float>(y: F, r: F, k1: F, k2: F, k3: F) -> DenormDouble<F> {
    let y = SemiDouble::new(y);

    // t = y * (k3 + r)
    let s = SemiDouble::new_qadd11(k3, r);
    let t = y * s;

    // t = y * (k2 + y * (k3 + r))
    let s = SemiDouble::new_qadd12(k2, t);
    let t = y * s;

    // y * (k1 + y * (k2 + y * (k3 + r)))
    let s = SemiDouble::new_qadd12(k1, t);
    y * s
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_tgamma<F: Float + FloatMath>() {
        use crate::tgamma;

        assert_is_nan!(tgamma(F::NAN));
        assert_is_nan!(tgamma(F::neg_infinity()));
        assert_is_nan!(tgamma(-F::one()));
        assert_is_nan!(tgamma(-F::two()));
        assert_is_nan!(tgamma(-F::largest()));
        assert_total_eq!(tgamma(F::INFINITY), F::INFINITY);
        assert_total_eq!(tgamma(F::ZERO), F::INFINITY);
        assert_total_eq!(tgamma(-F::ZERO), F::neg_infinity());
        assert_total_eq!(tgamma(F::one()), F::one());
        assert_total_eq!(tgamma(F::two()), F::one());
    }

    fn test_lgamma<F: Float + FloatMath>() {
        use crate::lgamma;

        let test_nan = |x: F| {
            let (r, sign) = lgamma(x);
            assert_is_nan!(r);
            assert_eq!(sign, 0);
        };
        let test_value = |x: F, r: F, sign: i8| {
            let (res, res_sign) = lgamma(x);
            assert_total_eq!(res, r);
            assert_eq!(res_sign, sign);
        };

        test_nan(F::NAN);
        test_nan(F::neg_infinity());
        test_value(F::INFINITY, F::INFINITY, 1);
        test_value(F::ZERO, F::INFINITY, 1);
        test_value(-F::ZERO, F::INFINITY, -1);
        test_value(-F::one(), F::INFINITY, 0);
        test_value(-F::two(), F::INFINITY, 0);
        test_value(-F::largest(), F::INFINITY, 0);
        test_value(F::one(), F::ZERO, 1);
        test_value(F::two(), F::ZERO, 1);
    }

    #[test]
    fn test_f32() {
        test_tgamma::<f32>();
        test_lgamma::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_tgamma::<crate::SoftF32>();
        test_lgamma::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_tgamma::<f64>();
        test_lgamma::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_tgamma::<crate::SoftF64>();
        test_lgamma::<crate::SoftF64>();
    }
}
