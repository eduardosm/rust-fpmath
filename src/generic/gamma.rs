use super::exp::{exp_split, hi_lo_exp_inner_common};
use super::ln::hi_lo_ln_hi_lo_inner;
use super::scalbn;
use super::sin_cos::{hi_lo_cos_inner, hi_lo_sin_inner};
use super::{Exp, Ln, ReduceHalfMulPi, SinCos, is_int, reduce_half_mul_pi};
use crate::double::{DenormDouble, NormDouble, SemiDouble};
use crate::traits::{Float, FloatConsts, Int as _};

pub(crate) trait Gamma: FloatConsts + SinCos + ReduceHalfMulPi + Exp + Ln {
    const LO_TH: Self;
    const HI_TH: Self;

    const TH_1: Self;
    const TH_2: Self;
    const TH_3: Self;

    const POLY_OFF: u8;

    const HALF_LN_2_PI: NormDouble<Self>;

    fn ln_gamma_poly_1(x: Self) -> (Self, Self, Self, Self);
    fn ln_gamma_poly_2(x: Self) -> (Self, Self, Self, Self);

    fn special_poly(x: Self) -> Self;
}

pub(crate) fn gamma<F: Gamma>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // gamma(±0) = ±inf
        F::INFINITY.copysign(x)
    } else if x >= F::HI_TH {
        // also handles x = inf
        F::INFINITY
    } else if e == F::MAX_RAW_EXP {
        // gamma(NaN or -inf) = NaN
        F::NAN
    } else if x.sign() && is_int(x) {
        // gamma(neg integer) = NaN
        F::NAN
    } else if x < F::LO_TH {
        // -inf and negative integers are handled above
        F::ZERO
    } else {
        gamma_inner(x)
    }
}

pub(crate) fn ln_gamma<F: Gamma>(x: F) -> (F, i8) {
    let e = x.raw_exp();
    let sign = x.sign();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // ln_gamma(0) = inf
        (F::INFINITY, if sign { -1 } else { 1 })
    } else if e == F::MAX_RAW_EXP {
        if !sign && x.raw_mant() == F::Raw::ZERO {
            // ln_gamma(inf) = inf
            (F::INFINITY, 1)
        } else {
            // ln_gamma(NaN or -inf) = NaN
            (F::NAN, 0)
        }
    } else if x.sign() && is_int(x) {
        // ln_gamma(neg integer) = inf
        (F::INFINITY, 0)
    } else if x == F::ONE || x == F::TWO {
        // ln_gamma(1 or 2) = 0
        // ensure positive zero
        (F::ZERO, 1)
    } else {
        ln_gamma_inner(x)
    }
}

fn gamma_inner<F: Gamma>(x: F) -> F {
    let (y, s) = gamma_inner_common(x);

    let (k, r_hi, r_lo) = exp_split(y.hi());
    let r_lo = r_lo + y.lo();
    let exp_y = hi_lo_exp_inner_common(r_hi, r_lo);

    scalbn((s.to_semi() * exp_y.to_semi()).to_single(), k)
}

fn ln_gamma_inner<F: Gamma>(x: F) -> (F, i8) {
    let (y, s) = gamma_inner_common(x);

    if y.hi() == F::INFINITY {
        (F::INFINITY, 1)
    } else {
        let s = s.to_norm();
        let (sign, abs_s) = if s.hi().sign() { (-1, -s) } else { (1, s) };

        let ln_s = hi_lo_ln_hi_lo_inner(abs_s, F::Exp::ZERO);

        (ln_s.qadd2(y).to_single(), sign)
    }
}

/// Returns `(y, s)` such as `Γ(x) = s * exp(y)`.
fn gamma_inner_common<F: Gamma>(x: F) -> (DenormDouble<F>, DenormDouble<F>) {
    // For x < 0.5, use gamma reflection formula:
    // Γ(x)*Γ(1-x) = π/sin(πx) => Γ(x) = π/(sin(πx)*Γ(1-x))
    let reflect = (x < F::HALF).then(|| {
        let (n, z) = reduce_half_mul_pi(x);
        let sinpix = match n {
            0 => hi_lo_sin_inner(z),
            1 => hi_lo_cos_inner(z),
            2 => -hi_lo_sin_inner(z),
            3 => -hi_lo_cos_inner(z),
            _ => unreachable!(),
        };
        // π / sin(πx)
        F::PI_EX / sinpix.to_semi()
    });
    // nx is always greater or equal to 0.5
    let nx = if reflect.is_some() {
        DenormDouble::new_sub11(F::ONE, x)
    } else {
        DenormDouble::new(x, F::ZERO)
    };

    // Based on the algorithm used in SLEEF.

    if nx.hi() < F::TH_2 {
        // For small values of `nx`, ln(Γ(nx)) is calculated using a polynomial.

        let nx = nx.hi();
        let y;
        let r;
        let k1;
        let k2;
        let k3;
        if nx < F::TH_1 {
            y = nx - F::ONE;
            (r, k1, k2, k3) = F::ln_gamma_poly_1(y);
        } else {
            y = nx - F::TWO;
            (r, k1, k2, k3) = F::ln_gamma_poly_2(y);
        };
        // r = ln(Γ(nx))
        let r = finish_poly(y, r, k1, k2, k3);

        if let Some(reflect) = reflect {
            // -ln(Γ(1 - x)), π / sin(πx)
            (-r, reflect)
        } else {
            // ln(Γ(x)), 1
            (r, DenormDouble::ONE)
        }
    } else {
        // For larger values of `nx`:
        // t = nx or nx + POLY_OFF
        // Γ(nx) = (P(1 / t) / t + 1) * t^(t - 0.5) * e^(-t) * √(2π)
        // P is a polynomial.

        let low = nx.hi() < F::TH_3;
        let t = if low {
            nx + F::cast_from(F::POLY_OFF)
        } else {
            nx
        };
        let tinv = F::ONE / t.to_single();

        // p = P(1 / t) * (1 / t) + 1
        let p1 = F::special_poly(tinv);
        let p = SemiDouble::new(p1) * SemiDouble::new(tinv) + F::ONE;

        // r = (t - 0.5) * ln(t) - t + 0.5 * ln(2π)
        //   = t * (ln(t) - 1) - 0.5 * ln(t) + 0.5 * ln(2π)
        let ln_t = hi_lo_ln_hi_lo_inner(t.to_norm(), F::Exp::ZERO);
        let r = t.to_semi() * (ln_t - F::ONE).to_semi() - ln_t.pmul1(F::HALF)
            + F::HALF_LN_2_PI.to_denorm();

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
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_gamma<F: Float + FloatMath>() {
        use crate::gamma;

        assert_is_nan!(gamma(F::NAN));
        assert_is_nan!(gamma(F::NEG_INFINITY));
        assert_is_nan!(gamma(-F::ONE));
        assert_is_nan!(gamma(-F::TWO));
        assert_is_nan!(gamma(-F::LARGEST));
        assert_total_eq!(gamma(F::INFINITY), F::INFINITY);
        assert_total_eq!(gamma(F::ZERO), F::INFINITY);
        assert_total_eq!(gamma(-F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(gamma(F::ONE), F::ONE);
        assert_total_eq!(gamma(F::TWO), F::ONE);
    }

    fn test_ln_gamma<F: Float + FloatMath>() {
        use crate::ln_gamma;

        let test_nan = |x: F| {
            let (r, sign) = ln_gamma(x);
            assert_is_nan!(r);
            assert_eq!(sign, 0);
        };
        let test_value = |x: F, r: F, sign: i8| {
            let (res, res_sign) = ln_gamma(x);
            assert_total_eq!(res, r);
            assert_eq!(res_sign, sign);
        };

        test_nan(F::NAN);
        test_nan(F::NEG_INFINITY);
        test_value(F::INFINITY, F::INFINITY, 1);
        test_value(F::ZERO, F::INFINITY, 1);
        test_value(-F::ZERO, F::INFINITY, -1);
        test_value(-F::ONE, F::INFINITY, 0);
        test_value(-F::TWO, F::INFINITY, 0);
        test_value(-F::LARGEST, F::INFINITY, 0);
        test_value(F::ONE, F::ZERO, 1);
        test_value(F::TWO, F::ZERO, 1);
    }

    #[test]
    fn test_f32() {
        test_gamma::<f32>();
        test_ln_gamma::<f32>();
    }

    #[test]
    fn test_f64() {
        test_gamma::<f64>();
        test_ln_gamma::<f64>();
    }
}
