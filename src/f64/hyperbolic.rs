use super::f64x2::F64x2;
use crate::generic::round_fi;
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0

// GENERATE: consts F64x2 LN_2
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1

impl crate::generic::Hyperbolic for f64 {
    fn sinh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp <= -600 {
            x
        } else if x_exp >= 10 {
            f64::INFINITY.copysign(x)
        } else {
            let (sinh, _, k) = sinh_cosh_inner(x);
            sinh.scalbn_to_f64(k)
        }
    }

    fn cosh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp <= -600 {
            1.0
        } else if x_exp >= 10 {
            f64::INFINITY
        } else {
            let (_, cosh, k) = sinh_cosh_inner(x);
            cosh.scalbn_to_f64(k)
        }
    }

    fn sinh_cosh_finite(x: Self) -> (Self, Self) {
        let x_exp = x.exponent();
        if x_exp <= -600 {
            (x, 1.0)
        } else if x_exp >= 10 {
            (f64::INFINITY.copysign(x), f64::INFINITY)
        } else {
            let (sinh, cosh, k) = sinh_cosh_inner(x);
            (sinh.scalbn_to_f64(k), cosh.scalbn_to_f64(k))
        }
    }

    fn tanh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp <= -600 {
            x
        } else if x_exp >= 8 {
            1.0f64.copysign(x)
        } else {
            let (sinh, cosh, _) = sinh_cosh_inner(x);
            (sinh / cosh).to_f64()
        }
    }
}

fn sinh_cosh_inner(x: f64) -> (F64x2, F64x2, i32) {
    #[inline]
    fn exp_m1_core(r: F64x2) -> (F64x2, F64x2) {
        // exp(r) = exp(r / s)^s
        // exp(-r) = exp(-r / s)^s
        // s = 2^5

        let rs = r.scalbn_fast(-5);

        // GENERATE: exp_m1_poly F64x2 6 -0.0108305 0.0108305
        const K2: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0x3C542E300EA57002); // 5.000000000000000043759312134915e-1
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555555, 0x3C6891CC6ECCB84F); // 1.666666666666666680701413615528e-1
        const K4: F64x2 = F64x2::from_bits(0x3FA555555554AD46, 0xBC4366ECBFFF15AE); // 4.166666666636813021769975331507e-2
        const K5: F64x2 = F64x2::from_bits(0x3F81111111108FA0, 0x3C274C442370496B); // 8.333333333275850317352901717945e-3
        const K6: F64x2 = F64x2::from_bits(0x3F56C1718E1A0E5A, 0x3BFDB83F67035EB1); // 1.388893979680802713641418637913e-3
        const K7: F64x2 = F64x2::from_bits(0x3F2A01A6362A32F7, 0xBBC5CBE5866A78EC); // 1.984134097096285320280180755654e-4

        let rs2 = rs.square();
        let mrs = -rs;
        let mut t1 = rs + horner!(rs2, rs, [K2, K3, K4, K5, K6, K7]);
        let mut t2 = mrs + horner!(rs2, mrs, [K2, K3, K4, K5, K6, K7]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..5 {
            t1 = t1.square() + t1.twice();
            t2 = t2.square() + t2.twice();
        }

        (t1, t2)
    }

    // x = k*ln(2) + r
    // k is an integer
    // |r| <= 0.5*ln(2)
    let (kf, k) = round_fi(x * LOG2_E);
    let r = x - kf * LN_2;

    // t1 = exp(r) - 1
    // t2 = exp(-r) - 1
    let (t1a, t1b) = exp_m1_core(r);

    // sinh(x) = (exp(x) - exp(-x)) / 2
    // cosh(x) = (exp(x) + exp(-x)) / 2

    if k == 0 {
        let s = (t1a - t1b).halve();
        let c = (t1a + t1b).halve() + 1.0;
        (s, c, 0)
    } else {
        let k = k as i32;
        let (k, ka, kb) = if k < 0 {
            (-k - 1, 2 * k, 0)
        } else {
            (k - 1, 0, -2 * k)
        };
        let t2a = (t1a + 1.0).scalbn_medium(ka);
        let t2b = (t1b + 1.0).scalbn_medium(kb);

        let s = t2a - t2b;
        let c = t2a + t2b;
        (s, c, k)
    }
}
