use crate::traits::Float as _;

// GENERATE: consts f64 FRAC_180_PI FRAC_1_PI
const FRAC_180_PI: f64 = f64::from_bits(0x404CA5DC1A63C1F8); // 5.729577951308232e1
const FRAC_1_PI: f64 = f64::from_bits(0x3FD45F306DC9C883); // 3.183098861837907e-1

impl crate::generic::InvTrigonometric for f32 {
    // GENERATE: consts f32 PI FRAC_PI_2
    const PI: f32 = f32::from_bits(0x40490FDB); // 3.1415927e0
    const FRAC_PI_2: f32 = f32::from_bits(0x3FC90FDB); // 1.5707964e0

    fn asin_finite(x: Self) -> Self {
        asin_core(x) as f32
    }

    fn acos_finite(x: Self) -> Self {
        acos_core(x) as f32
    }

    fn atan_finite(x: Self) -> Self {
        atan_core(x) as f32
    }

    fn atan2_finite(y: Self, x: Self) -> Self {
        atan2_core(y, x) as f32
    }

    fn asind_finite(x: Self) -> Self {
        (asin_core(x) * FRAC_180_PI) as f32
    }

    fn acosd_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_180_PI) as f32
    }

    fn atand_finite(x: Self) -> Self {
        (atan_core(x) * FRAC_180_PI) as f32
    }

    fn atan2d_finite(y: Self, x: Self) -> Self {
        (atan2_core(y, x) * FRAC_180_PI) as f32
    }

    fn asinpi_finite(x: Self) -> Self {
        (asin_core(x) * FRAC_1_PI) as f32
    }

    fn acospi_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_1_PI) as f32
    }

    fn atanpi_finite(x: Self) -> Self {
        (atan_core(x) * FRAC_1_PI) as f32
    }

    fn atan2pi_finite(y: Self, x: Self) -> Self {
        (atan2_core(y, x) * FRAC_1_PI) as f32
    }
}

fn asin_core(x: f32) -> f64 {
    // GENERATE: consts f64 FRAC_PI_2
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0

    #[inline]
    fn asin_poly(x: f64, x2: f64, x3: f64) -> f64 {
        // GENERATE: asin_poly f64 7
        const K3: f64 = f64::from_bits(0x3FC555556DCEC097); // 1.6666667806339738e-1
        const K5: f64 = f64::from_bits(0x3FB33324196CC1B0); // 7.499909992907905e-2
        const K7: f64 = f64::from_bits(0x3FA6DECF11049302); // 4.466864664781235e-2
        const K9: f64 = f64::from_bits(0x3F9EBD85687015AF); // 3.0019840716723752e-2
        const K11: f64 = f64::from_bits(0x3F99B8A837AD73CC); // 2.511847343281541e-2
        const K13: f64 = f64::from_bits(0x3F78EBA7BEDCFA90); // 6.0841133652603935e-3
        const K15: f64 = f64::from_bits(0x3FA28A536EDECDA2); // 3.621159294507527e-2

        x + horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15])
    }

    let x = f64::from(x);
    let xexp = x.exponent();
    if xexp < -100 {
        // zero or tiny
        x
    } else if xexp < -1 {
        // |x| < 0.5
        let x2 = x * x;
        let x3 = x2 * x;

        asin_poly(x, x2, x3)
    } else {
        // |x| >= 0.5
        // |asin(x)| = π/2 - 2 * asin(sqrt((1 - |x|) / 2))

        // y = sqrt((1 - |x|) / 2)
        let y2 = (1.0 - x.abs()) * 0.5;
        let y = crate::f64::fast_sqrt(y2).0;
        let y3 = y2 * y;

        // t1 = asin(y)
        let t1 = asin_poly(y, y2, y3);

        // t3 = π/2 - 2 * asin(sqrt((1 - |x|) / 2)) = π/2 - 2 * t1
        let t3 = FRAC_PI_2 - 2.0 * t1;

        t3.copysign(x)
    }
}

fn acos_core(x: f32) -> f64 {
    // GENERATE: consts f64 FRAC_PI_2 PI
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
    const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0

    #[inline]
    fn asin_poly(x: f64, x2: f64, x3: f64) -> f64 {
        // GENERATE: asin_poly f64 7
        const K3: f64 = f64::from_bits(0x3FC555556DCEC097); // 1.6666667806339738e-1
        const K5: f64 = f64::from_bits(0x3FB33324196CC1B0); // 7.499909992907905e-2
        const K7: f64 = f64::from_bits(0x3FA6DECF11049302); // 4.466864664781235e-2
        const K9: f64 = f64::from_bits(0x3F9EBD85687015AF); // 3.0019840716723752e-2
        const K11: f64 = f64::from_bits(0x3F99B8A837AD73CC); // 2.511847343281541e-2
        const K13: f64 = f64::from_bits(0x3F78EBA7BEDCFA90); // 6.0841133652603935e-3
        const K15: f64 = f64::from_bits(0x3FA28A536EDECDA2); // 3.621159294507527e-2

        x + horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15])
    }

    let x = f64::from(x);

    // acos(x) = π/2 - asin(x)
    if x.exponent() < -1 {
        // |x| < 0.5
        let x2 = x * x;
        let x3 = x2 * x;

        // t1 = asin(x)
        let t1 = asin_poly(x, x2, x3);

        // acos(x) = π/2 - asin(x) = π/2 - t1
        FRAC_PI_2 - t1
    } else {
        // |x| >= 0.5
        // |asin(x)| = π/2 - 2 * asin(sqrt((1 - |x|) / 2))

        // y = sqrt((1 - |x|) / 2)
        let y2 = (1.0 - x.abs()) * 0.5;
        let y = crate::f64::fast_sqrt(y2).0;
        let y3 = y2 * y;

        // t1 = asin(y)
        let t1 = asin_poly(y, y2, y3);

        // t2 = 2 * asin(y) = 2 * t1
        let t2 = 2.0 * t1;

        if x.is_sign_negative() {
            // acos(x) = π/2 + |asin(x)|
            //         = π/2 + (π/2 - 2 * asin(y))
            //         = π - 2 * asin(y)
            //         = π - t2
            PI - t2
        } else {
            // acos(x) = π/2 - |asin(x)|
            //         = π/2 - (π/2 - 2 * asin(y))
            //         = 2 * asin(y)
            //         = t2
            t2
        }
    }
}

fn atan_core(x: f32) -> f64 {
    // GENERATE: consts f64 FRAC_PI_2
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0

    #[inline]
    fn atan_poly(x: f64, x2: f64, x3: f64) -> f64 {
        // GENERATE: atan_poly f64 12 -0.00001 1.0
        const K3: f64 = f64::from_bits(0xBFD555554F68C559); // -3.333333278161113e-1
        const K5: f64 = f64::from_bits(0x3FC99997094AB706); // 1.9999969438270443e-1
        const K7: f64 = f64::from_bits(0xBFC248EE9E86D65D); // -1.428507112556642e-1
        const K9: f64 = f64::from_bits(0x3FBC6D143BCF5365); // 1.1103941402627766e-1
        const K11: f64 = f64::from_bits(0xBFB725C6FE9111C5); // -9.042018618590138e-2
        const K13: f64 = f64::from_bits(0x3FB3205A97388D0B); // 7.471243087688977e-2
        const K15: f64 = f64::from_bits(0xBFAE8F4BC9F39CA8); // -5.96870121024094e-2
        const K17: f64 = f64::from_bits(0x3FA5EFE28EEACC68); // 4.284580225557805e-2
        const K19: f64 = f64::from_bits(0xBF99D06CD85E3CD4); // -2.5209141450940845e-2
        const K21: f64 = f64::from_bits(0x3F865E4AB262F61C); // 1.0922034806180973e-2
        const K23: f64 = f64::from_bits(0xBF68A9059ACF59C6); // -3.0102834193528076e-3
        const K25: f64 = f64::from_bits(0x3F3985DFF09C577A); // 3.894492843853275e-4

        x + horner!(
            x3,
            x2,
            [K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25]
        )
    }

    let x = f64::from(x);
    if x.exponent() < -100 {
        // zero or tiny
        x
    } else if x.abs() <= 1.0 {
        let x2 = x * x;
        let x3 = x2 * x;
        atan_poly(x, x2, x3)
    } else {
        // atan(x) = ±pi/2 - atan(1 / x)
        let inv_x = x.recip();
        let inv_x2 = inv_x * inv_x;
        let inv_x3 = inv_x2 * inv_x;

        FRAC_PI_2.copysign(x) - atan_poly(inv_x, inv_x2, inv_x3)
    }
}

fn atan2_core(n: f32, d: f32) -> f64 {
    // GENERATE: consts f64 FRAC_PI_2
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0

    #[inline]
    fn atan_poly(x: f64, x2: f64, x3: f64) -> f64 {
        // GENERATE: atan_poly f64 12 -0.00001 1.0
        const K3: f64 = f64::from_bits(0xBFD555554F68C559); // -3.333333278161113e-1
        const K5: f64 = f64::from_bits(0x3FC99997094AB706); // 1.9999969438270443e-1
        const K7: f64 = f64::from_bits(0xBFC248EE9E86D65D); // -1.428507112556642e-1
        const K9: f64 = f64::from_bits(0x3FBC6D143BCF5365); // 1.1103941402627766e-1
        const K11: f64 = f64::from_bits(0xBFB725C6FE9111C5); // -9.042018618590138e-2
        const K13: f64 = f64::from_bits(0x3FB3205A97388D0B); // 7.471243087688977e-2
        const K15: f64 = f64::from_bits(0xBFAE8F4BC9F39CA8); // -5.96870121024094e-2
        const K17: f64 = f64::from_bits(0x3FA5EFE28EEACC68); // 4.284580225557805e-2
        const K19: f64 = f64::from_bits(0xBF99D06CD85E3CD4); // -2.5209141450940845e-2
        const K21: f64 = f64::from_bits(0x3F865E4AB262F61C); // 1.0922034806180973e-2
        const K23: f64 = f64::from_bits(0xBF68A9059ACF59C6); // -3.0102834193528076e-3
        const K25: f64 = f64::from_bits(0x3F3985DFF09C577A); // 3.894492843853275e-4

        x + horner!(
            x3,
            x2,
            [K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25]
        )
    }

    let mut n = f64::from(n);
    let mut d = f64::from(d);

    let ysgn = n.is_sign_negative();
    let xsgn = d.is_sign_negative();

    let mut off = 0.0;
    if xsgn {
        off = 2.0f64.set_sign(ysgn);
    }
    if n.abs() > d.abs() {
        (n, d) = (d, n);
        n = -n;
        off += 1.0f64.set_sign(ysgn ^ xsgn);
    }

    let z = n / d;
    if off == 0.0 && z.exponent() <= -32 {
        // atan2(y, x) ~= y/x = n/d
        z
    } else {
        // atan2(y, x) = atan(n/d) + off * π/2
        let z2 = z * z;
        let z3 = z2 * z;
        atan_poly(z, z2, z3) + off * FRAC_PI_2
    }
}
