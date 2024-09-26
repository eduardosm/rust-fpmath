use super::exp::exp_split;
use super::{scalbn_medium, Exp};
use crate::double::DenormDouble;
use crate::traits::{CastInto as _, Float, Int as _, Like};

pub(crate) trait SinhCosh<L = Like<Self>>: Exp {
    fn expo2_hi_th() -> Self;
}

pub(crate) fn sinh<F: SinhCosh>(x: F) -> F {
    let e = x.raw_exp();
    if x >= F::expo2_hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if x <= -F::expo2_hi_th() {
        // also handles x = -inf
        F::neg_infinity()
    } else if e == F::MAX_RAW_EXP || e <= F::RawExp::ONE {
        // propagate NaN
        // or
        // very small, includes subnormal and zero
        // sinh(x) ~= x
        // also handles sinh(-0) = -0
        x
    } else {
        let (s, _) = sinh_cosh_inner(x);
        s
    }
}

pub(crate) fn cosh<F: SinhCosh>(x: F) -> F {
    let e = x.raw_exp();
    if x.abs() >= F::expo2_hi_th() {
        // also handles x = ±inf
        F::INFINITY
    } else if e == F::MAX_RAW_EXP {
        // Propagate NaN
        x
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // cosh(x) ~= 1
        F::one()
    } else {
        let (_, c) = sinh_cosh_inner(x);
        c
    }
}

pub(crate) fn sinh_cosh<F: SinhCosh>(x: F) -> (F, F) {
    let e = x.raw_exp();
    if x >= F::expo2_hi_th() {
        // also handles x = inf
        (F::INFINITY, F::INFINITY)
    } else if x <= -F::expo2_hi_th() {
        // also handles x = -inf
        (F::neg_infinity(), F::INFINITY)
    } else if e == F::MAX_RAW_EXP {
        // Propagate NaN
        (x, x)
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // sinh(x) ~= x
        // cosh(x) ~= 1
        // also handles sinh(-0) = -0
        (x, F::one())
    } else {
        let (s, c) = sinh_cosh_inner(x);
        (s, c)
    }
}

fn sinh_cosh_inner<F: Exp>(x: F) -> (F, F) {
    // Split |x| into k, r_hi, r_lo such as:
    //  - |x| = k*ln(2) + r_hi + r_lo
    //  - k is an integer
    //  - |r_hi| <= 0.5*ln(2)
    let absx = x.abs();
    let (k, r_hi, r_lo) = exp_split(absx);
    let (r_hi, r_lo) = F::norm_hi_lo_full(r_hi, r_lo);

    // t1a = exp(r_hi + r_lo) - 1 - r_hi
    // t1b = exp(-r_hi - r_lo) - 1 + r_hi
    let (t1a, t1b) = sinh_cosh_inner_common_1(r_hi, r_lo);

    if k > F::MANT_BITS.into() {
        let t2 = scalbn_medium((r_hi + t1a) + F::one(), k - 1);
        (t2.copysign(x), t2)
    } else {
        // abss = |sinh(x)| = (exp(|x|) - exp(-|x|)) / 2
        // c = cosh(x) = (exp(|x|) + exp(-|x|)) / 2
        let (abss, c) = sinh_cosh_inner_common_2(k, r_hi, t1a, t1b);

        (abss.to_single().copysign(x), c.to_single())
    }
}

/// Returns `(ta, tb)`
///
/// Such as
/// * `ta = exp(r_hi + r_lo) - 1 - r_hi`
/// * `tb = exp(-r_hi - r_lo) - 1 + r_hi`
pub(super) fn sinh_cosh_inner_common_1<F: Exp>(r_hi: F, r_lo: F) -> (F, F) {
    // pseudo-consts
    let three = F::one() + F::two();
    let six = three * F::two();

    // Similar as done in `exp_m1_inner` in exp.rs, but to calculate
    // both exp(r) and exp(-r)
    let r2 = r_hi * r_hi;
    let hr = F::half() * r_hi;
    let hr2 = F::half() * r2;

    let t1 = F::exp_m1_special_poly(r2);

    let t2a = three - t1 * hr;
    let t3a = hr2 * ((t1 - t2a) / (six - r_hi * t2a));

    let t2b = three + t1 * hr;
    let t3b = hr2 * ((t1 - t2b) / (six + r_hi * t2b));

    // t4a = exp(r_hi + r_lo) - 1 - r
    let t4a = (r_hi * (r_lo - t3a) + r_lo) + hr2;

    // t4b = exp(-r_hi - r_lo) - 1 + r
    let t4b = (r_hi * (r_lo + t3b) - r_lo) + hr2;

    (t4a, t4b)
}

/// Calculates `(abss, c)`
///
/// Such as
/// * `abs = (1 + r + t1a) * 2^(k - 1) - (1 - r + t1b) * 2^(-k - 1)`
/// * `c = (1 + r + t1a) * 2^(k - 1) + (1 - r + t1b) * 2^(-k - 1)`
pub(super) fn sinh_cosh_inner_common_2<F: Float>(
    k: i32,
    r: F,
    t1a: F,
    t1b: F,
) -> (DenormDouble<F>, DenormDouble<F>) {
    let s1a = F::exp2i_fast((k - 1).cast_into());
    let sra = r * s1a;
    let st1a = t1a * s1a;

    let s1b = F::exp2i_fast((-k - 1).cast_into());
    let srb = r * s1b;
    let st1b = t1b * s1b;

    if k <= 1 {
        let t2a = s1a - s1b;
        let t3a = DenormDouble::new_qadd11(sra, srb).qradd1(t2a);
        let abss = t3a.ladd(st1a - st1b);

        let t2b = s1a + s1b;
        let t3b = DenormDouble::new_qsub11(sra, srb).qradd1(t2b);
        let c = t3b.ladd(st1a + st1b);

        (abss, c)
    } else {
        let t2 = DenormDouble::new_qadd11(sra, st1a).qradd1(s1a);

        let t3a = t2.qsub1(s1b);
        let abss = t3a.ladd(srb - st1b);

        let t3b = t2.qadd1(s1b);
        let c = t3b.ladd(st1b - srb);

        (abss, c)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>(hi_th: &str) {
        use crate::{cosh, sinh, sinh_cosh};

        let test_nan = |arg: F| {
            let sinh1 = sinh(arg);
            let cosh1 = cosh(arg);
            let (sinh2, cosh2) = sinh_cosh(arg);
            assert_is_nan!(sinh1);
            assert_is_nan!(cosh1);
            assert_is_nan!(sinh2);
            assert_is_nan!(cosh2);
        };

        let test_value = |arg: F, expected_sinh: F, expected_cosh: F| {
            let sinh1 = sinh(arg);
            let cosh1 = cosh(arg);
            let (sinh2, cosh2) = sinh_cosh(arg);
            assert_total_eq!(sinh1, expected_sinh);
            assert_total_eq!(cosh1, expected_cosh);
            assert_total_eq!(sinh2, expected_sinh);
            assert_total_eq!(cosh2, expected_cosh);
        };

        let hi_th = F::parse(hi_th);

        test_nan(F::NAN);
        test_value(F::INFINITY, F::INFINITY, F::INFINITY);
        test_value(F::neg_infinity(), F::neg_infinity(), F::INFINITY);
        test_value(hi_th, F::INFINITY, F::INFINITY);
        test_value(-hi_th, F::neg_infinity(), F::INFINITY);
        test_value(hi_th + F::half(), F::INFINITY, F::INFINITY);
        test_value(-(hi_th + F::half()), F::neg_infinity(), F::INFINITY);
        test_value(hi_th + F::one(), F::INFINITY, F::INFINITY);
        test_value(-(hi_th + F::one()), F::neg_infinity(), F::INFINITY);
        test_value(F::ZERO, F::ZERO, F::one());
        test_value(-F::ZERO, -F::ZERO, F::one());
    }

    #[test]
    fn test_f32() {
        test::<f32>("89.5");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>("89.5");
    }

    #[test]
    fn test_f64() {
        test::<f64>("710.5");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>("710.5");
    }
}
