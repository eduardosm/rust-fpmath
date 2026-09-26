use super::f64x2::F64x2;
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_180_PI FRAC_1_PI
const FRAC_180_PI: F64x2 = F64x2::from_bits(0x404CA5DC1A63C1F8, 0xBCE1E7AB456405F9); // 5.729577951308232087679815481411e1
const FRAC_1_PI: F64x2 = F64x2::from_bits(0x3FD45F306DC9C883, 0xBC76B01EC5417056); // 3.183098861837906715377675267450e-1

impl crate::generic::InvTrigonometric for f64 {
    // GENERATE: consts f64 PI FRAC_PI_2
    const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0

    fn asin_finite(x: Self) -> Self {
        asin_core(x).to_f64()
    }

    fn acos_finite(x: Self) -> Self {
        acos_core(x).to_f64()
    }

    fn atan_finite(x: Self) -> Self {
        atan_core(x).to_f64()
    }

    fn atan2_finite(y: Self, x: Self) -> Self {
        let (r, edelta) = atan2_core(y, x);
        r.scalbn_to_f64(edelta)
    }

    fn asind_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_180_PI)
                .scalbn_to_f64(-106)
                .copysign(x)
        } else {
            (asin_core(x) * FRAC_180_PI).to_f64()
        }
    }

    fn acosd_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_180_PI).to_f64()
    }

    fn atand_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_180_PI)
                .scalbn_to_f64(-106)
                .copysign(x)
        } else {
            (atan_core(x) * FRAC_180_PI).to_f64()
        }
    }

    fn atan2d_finite(y: Self, x: Self) -> Self {
        // `F64x2` can lose the sign of zero, `copysign` recovers it.
        let (r, edelta) = atan2_core(y, x);
        (r * FRAC_180_PI).scalbn_to_f64(edelta).copysign(y)
    }

    fn asinpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_1_PI)
                .scalbn_to_f64(-106)
                .copysign(x)
        } else {
            (asin_core(x) * FRAC_1_PI).to_f64()
        }
    }

    fn acospi_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_1_PI).to_f64()
    }

    fn atanpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_1_PI)
                .scalbn_to_f64(-106)
                .copysign(x)
        } else {
            (atan_core(x) * FRAC_1_PI).to_f64()
        }
    }

    fn atan2pi_finite(y: Self, x: Self) -> Self {
        // `F64x2` can lose the sign of zero, `copysign` recovers it.
        let (r, edelta) = atan2_core(y, x);
        (r * FRAC_1_PI).scalbn_to_f64(edelta).copysign(y)
    }
}

fn asin_core(x: f64) -> F64x2 {
    // GENERATE: consts F64x2 FRAC_PI_2
    const FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FF921FB54442D18, 0x3C91A62633145C07); // 1.570796326794896619231321691640e0

    #[inline]
    fn asin_poly(x: F64x2, x2: F64x2, x3: F64x2) -> F64x2 {
        // GENERATE: asin_poly F64x2 15
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555556, 0x3C3FD402B9ACC12D); // 1.666666666666666868957921116186e-1
        const K5: F64x2 = F64x2::from_bits(0x3FB3333333333171, 0xBC48AA48BF2C58F6); // 7.499999999999374954572555801348e-2
        const K7: F64x2 = F64x2::from_bits(0x3FA6DB6DB6DD0B75, 0x3C490DC78356F66C); // 4.464285714359210032926867267681e-2
        const K9: F64x2 = F64x2::from_bits(0x3F9F1C71C652A920, 0xBC3C9999C6564CCE); // 3.038194439856411866287044220888e-2
        const K11: F64x2 = F64x2::from_bits(0x3F96E8BA4CB5D628, 0xBC319E3998AA4880); // 2.237216084673412834742706631679e-2
        const K13: F64x2 = F64x2::from_bits(0x3F91C4E95094DC2A, 0xBC2639DAD7993A8D); // 1.735271982508038791356636690954e-2
        const K15: F64x2 = F64x2::from_bits(0x3F8C9A037B7E7D9A, 0x3C2656B3CC2DA22E); // 1.396563263534638240679125901867e-2
        const K17: F64x2 = F64x2::from_bits(0x3F87A3392B6C541C, 0x3C2DBC88A34A3A1A); // 1.154179252402270314358000559656e-2
        const K19: F64x2 = F64x2::from_bits(0x3F842E6D8185255E, 0x3C048B5CB967E668); // 9.854178919925512536873675592482e-3
        const K21: F64x2 = F64x2::from_bits(0x3F7FCBCA1E90290B, 0xBBFF985133F1BDA5); // 7.762708214622780486335590695343e-3
        const K23: F64x2 = F64x2::from_bits(0x3F855455B99746F3, 0xBC14CF4EB0C745C2); // 1.041476223591384301755867599605e-2
        const K25: F64x2 = F64x2::from_bits(0xBF72823A4CCD4799, 0xBC10DF34B9879994); // -4.518726095584479615314454218699e-3
        const K27: F64x2 = F64x2::from_bits(0x3FA0918BE46E6C7B, 0xBC273AAD7AF7555A); // 3.236043132275142418006266820483e-2
        const K29: F64x2 = F64x2::from_bits(0xBFA2462C2F2A1216, 0x3C28CDC3A93B0C4D); // -3.569162441002419960338541647919e-2
        const K31: F64x2 = F64x2::from_bits(0x3FA2EDBCCAD8E19B, 0xBC1675F3900C6760); // 3.697004295503632906407903775709e-2

        x + horner!(
            x3,
            x2,
            [
                K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31
            ]
        )
    }

    let x = F64x2::new1(x);
    if x.hi().exponent() < -1 {
        // |x| < 0.5
        let x2 = x * x;
        let x3 = x2 * x;

        asin_poly(x, x2, x3)
    } else {
        // |x| >= 0.5
        // |asin(x)| = π/2 - 2 * asin(sqrt((1 - |x|) / 2))

        // y = sqrt((1 - |x|) / 2)
        let y2 = (1.0 - x.abs()).halve();
        let y = y2.sqrt();
        let y3 = y2 * y;

        // t1 = asin(y)
        let t1 = asin_poly(y, y2, y3);

        // t3 = π/2 - 2 * asin(sqrt((1 - |x|) / 2)) = π/2 - 2 * t1
        let t3 = FRAC_PI_2 - t1.twice();

        if x.hi().is_sign_negative() { -t3 } else { t3 }
    }
}

fn acos_core(x: f64) -> F64x2 {
    // GENERATE: consts F64x2 FRAC_PI_2 PI
    const FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FF921FB54442D18, 0x3C91A62633145C07); // 1.570796326794896619231321691640e0
    const PI: F64x2 = F64x2::from_bits(0x400921FB54442D18, 0x3CA1A62633145C07); // 3.141592653589793238462643383280e0

    #[inline]
    fn asin_poly(x: F64x2, x2: F64x2, x3: F64x2) -> F64x2 {
        // GENERATE: asin_poly F64x2 15
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555556, 0x3C3FD402B9ACC12D); // 1.666666666666666868957921116186e-1
        const K5: F64x2 = F64x2::from_bits(0x3FB3333333333171, 0xBC48AA48BF2C58F6); // 7.499999999999374954572555801348e-2
        const K7: F64x2 = F64x2::from_bits(0x3FA6DB6DB6DD0B75, 0x3C490DC78356F66C); // 4.464285714359210032926867267681e-2
        const K9: F64x2 = F64x2::from_bits(0x3F9F1C71C652A920, 0xBC3C9999C6564CCE); // 3.038194439856411866287044220888e-2
        const K11: F64x2 = F64x2::from_bits(0x3F96E8BA4CB5D628, 0xBC319E3998AA4880); // 2.237216084673412834742706631679e-2
        const K13: F64x2 = F64x2::from_bits(0x3F91C4E95094DC2A, 0xBC2639DAD7993A8D); // 1.735271982508038791356636690954e-2
        const K15: F64x2 = F64x2::from_bits(0x3F8C9A037B7E7D9A, 0x3C2656B3CC2DA22E); // 1.396563263534638240679125901867e-2
        const K17: F64x2 = F64x2::from_bits(0x3F87A3392B6C541C, 0x3C2DBC88A34A3A1A); // 1.154179252402270314358000559656e-2
        const K19: F64x2 = F64x2::from_bits(0x3F842E6D8185255E, 0x3C048B5CB967E668); // 9.854178919925512536873675592482e-3
        const K21: F64x2 = F64x2::from_bits(0x3F7FCBCA1E90290B, 0xBBFF985133F1BDA5); // 7.762708214622780486335590695343e-3
        const K23: F64x2 = F64x2::from_bits(0x3F855455B99746F3, 0xBC14CF4EB0C745C2); // 1.041476223591384301755867599605e-2
        const K25: F64x2 = F64x2::from_bits(0xBF72823A4CCD4799, 0xBC10DF34B9879994); // -4.518726095584479615314454218699e-3
        const K27: F64x2 = F64x2::from_bits(0x3FA0918BE46E6C7B, 0xBC273AAD7AF7555A); // 3.236043132275142418006266820483e-2
        const K29: F64x2 = F64x2::from_bits(0xBFA2462C2F2A1216, 0x3C28CDC3A93B0C4D); // -3.569162441002419960338541647919e-2
        const K31: F64x2 = F64x2::from_bits(0x3FA2EDBCCAD8E19B, 0xBC1675F3900C6760); // 3.697004295503632906407903775709e-2

        x + horner!(
            x3,
            x2,
            [
                K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31
            ]
        )
    }

    let x = F64x2::new1(x);

    // acos(x) = π/2 - asin(x)
    if x.hi().exponent() < -1 {
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
        let y = y2.sqrt();
        let y3 = y2 * y;

        // t1 = asin(y)
        let t1 = asin_poly(y, y2, y3);

        // t2 = 2 * asin(y) = 2 * t1
        let t2 = t1.twice();

        if x.hi().is_sign_negative() {
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

fn atan_core(x: f64) -> F64x2 {
    // GENERATE: consts F64x2 FRAC_PI_2
    const FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FF921FB54442D18, 0x3C91A62633145C07); // 1.570796326794896619231321691640e0

    #[inline]
    fn atan_poly(x: F64x2, x2: F64x2, x3: F64x2) -> F64x2 {
        // GENERATE: atan_poly F64x2 23 -0.00001 1.0
        const K3: F64x2 = F64x2::from_bits(0xBFD5555555555554, 0xBC6BE0BF29E6CD8E); // -3.333333333333332714085842633092e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC9999999999809, 0x3C5C3CB33EBA06FD); // 1.999999999999888872393555938531e-1
        const K7: F64x2 = F64x2::from_bits(0xBFC249249248B8B9, 0xBC521AAB9492E989); // -1.428571428563765494512391719635e-1
        const K9: F64x2 = F64x2::from_bits(0x3FBC71C71C52771D, 0xBC5163647C4A78D7); // 1.111111110826325622012522842176e-1
        const K11: F64x2 = F64x2::from_bits(0xBFB745D17184A47F, 0xBC451028DBC118E4); // -9.090909024657099996796288877484e-2
        const K13: F64x2 = F64x2::from_bits(0x3FB3B13AE6867C0C, 0xBC56C17C99A06676); // 7.692306640682050117908280328992e-2
        const K15: F64x2 = F64x2::from_bits(0xBFB1110F0BAF08D5, 0xBC406ACB60837A87); // -6.666654620406416332038950930503e-2
        const K17: F64x2 = F64x2::from_bits(0x3FAE1DFB5E4EB8FE, 0x3C49C52656C0761F); // 5.882249380316472155595710125919e-2
        const K19: F64x2 = F64x2::from_bits(0xBFAAF1A000881FCF, 0x3C4A358D15B7870B); // -5.262470251551542154201320989941e-2
        const K21: F64x2 = F64x2::from_bits(0x3FA85CCCC06574D9, 0xBC35DBD62A4BE34C); // 4.758300636852092500315033964697e-2
        const K23: F64x2 = F64x2::from_bits(0xBFA62EE7B3754D6A, 0xBC390C3F5D5A1686); // -4.332660737659462104946945094290e-2
        const K25: F64x2 = F64x2::from_bits(0x3FA436D182F9E7D3, 0x3C06DADE4E63892B); // 3.948073123770225888395432340211e-2
        const K27: F64x2 = F64x2::from_bits(0xBFA236BDBC6C29E3, 0xBC31885C38546FF3); // -3.557389187888639848332458869199e-2
        const K29: F64x2 = F64x2::from_bits(0x3F9FCD44DEE81CB1, 0x3C2B4DD650D9BCFD); // 3.105647669329120596115913888384e-2
        const K31: F64x2 = F64x2::from_bits(0xBF9A237C144BAE1B, 0xBC3B0FEF76CA2078); // -2.552598832964579422410845993396e-2
        const K33: F64x2 = F64x2::from_bits(0x3F938B977542E09E, 0x3C3398B2898C14B1); // 1.908718732076463665199194660137e-2
        const K35: F64x2 = F64x2::from_bits(0xBF89A2F265CB6A84, 0xBC12A235126683DC); // -1.251782773163268744627397696142e-2
        const K37: F64x2 = F64x2::from_bits(0x3F7C65EC9711289A, 0xBC147B932323683D); // 6.933139972782020946059224346584e-3
        const K39: F64x2 = F64x2::from_bits(0xBF6980ABBC25DFED, 0xBC0374FF51B89672); // -3.113112850027014096722356271323e-3
        const K41: F64x2 = F64x2::from_bits(0x3F51ABA0D832F840, 0xBBDE640F1D3AD41C); // 1.078517031416051349183701491189e-3
        const K43: F64x2 = F64x2::from_bits(0xBF319D51BE7DBE7A, 0x3BD4A9E83110A63C); // -2.687763758250982140852720885938e-4
        const K45: F64x2 = F64x2::from_bits(0x3F066187B926004B, 0x3B57E8B70AE73E4B); // 4.268832630582914715025741411327e-5
        const K47: F64x2 = F64x2::from_bits(0xBECB236D2148CA77, 0x3B54350CDAE4FFCD); // -3.235147469271306583022107945664e-6

        x + horner!(
            x3,
            x2,
            [
                K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31, K33, K35,
                K37, K39, K41, K43, K45, K47
            ]
        )
    }

    let x = F64x2::new1(x);
    if x.hi().abs() <= 1.0 {
        let x2 = x * x;
        let x3 = x2 * x;
        atan_poly(x, x2, x3)
    } else {
        // atan(x) = ±pi/2 - atan(1 / x)
        let inv_x = x.recip();
        let inv_x2 = inv_x * inv_x;
        let inv_x3 = inv_x2 * inv_x;

        let t1 = atan_poly(inv_x, inv_x2, inv_x3);
        if x.hi().is_sign_negative() {
            -FRAC_PI_2 - t1
        } else {
            FRAC_PI_2 - t1
        }
    }
}

fn atan2_core(mut n: f64, mut d: f64) -> (F64x2, i32) {
    // GENERATE: consts F64x2 FRAC_PI_2
    const FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FF921FB54442D18, 0x3C91A62633145C07); // 1.570796326794896619231321691640e0

    #[inline]
    fn atan_poly(x: F64x2, x2: F64x2, x3: F64x2) -> F64x2 {
        // GENERATE: atan_poly F64x2 23 -0.00001 1.0
        const K3: F64x2 = F64x2::from_bits(0xBFD5555555555554, 0xBC6BE0BF29E6CD8E); // -3.333333333333332714085842633092e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC9999999999809, 0x3C5C3CB33EBA06FD); // 1.999999999999888872393555938531e-1
        const K7: F64x2 = F64x2::from_bits(0xBFC249249248B8B9, 0xBC521AAB9492E989); // -1.428571428563765494512391719635e-1
        const K9: F64x2 = F64x2::from_bits(0x3FBC71C71C52771D, 0xBC5163647C4A78D7); // 1.111111110826325622012522842176e-1
        const K11: F64x2 = F64x2::from_bits(0xBFB745D17184A47F, 0xBC451028DBC118E4); // -9.090909024657099996796288877484e-2
        const K13: F64x2 = F64x2::from_bits(0x3FB3B13AE6867C0C, 0xBC56C17C99A06676); // 7.692306640682050117908280328992e-2
        const K15: F64x2 = F64x2::from_bits(0xBFB1110F0BAF08D5, 0xBC406ACB60837A87); // -6.666654620406416332038950930503e-2
        const K17: F64x2 = F64x2::from_bits(0x3FAE1DFB5E4EB8FE, 0x3C49C52656C0761F); // 5.882249380316472155595710125919e-2
        const K19: F64x2 = F64x2::from_bits(0xBFAAF1A000881FCF, 0x3C4A358D15B7870B); // -5.262470251551542154201320989941e-2
        const K21: F64x2 = F64x2::from_bits(0x3FA85CCCC06574D9, 0xBC35DBD62A4BE34C); // 4.758300636852092500315033964697e-2
        const K23: F64x2 = F64x2::from_bits(0xBFA62EE7B3754D6A, 0xBC390C3F5D5A1686); // -4.332660737659462104946945094290e-2
        const K25: F64x2 = F64x2::from_bits(0x3FA436D182F9E7D3, 0x3C06DADE4E63892B); // 3.948073123770225888395432340211e-2
        const K27: F64x2 = F64x2::from_bits(0xBFA236BDBC6C29E3, 0xBC31885C38546FF3); // -3.557389187888639848332458869199e-2
        const K29: F64x2 = F64x2::from_bits(0x3F9FCD44DEE81CB1, 0x3C2B4DD650D9BCFD); // 3.105647669329120596115913888384e-2
        const K31: F64x2 = F64x2::from_bits(0xBF9A237C144BAE1B, 0xBC3B0FEF76CA2078); // -2.552598832964579422410845993396e-2
        const K33: F64x2 = F64x2::from_bits(0x3F938B977542E09E, 0x3C3398B2898C14B1); // 1.908718732076463665199194660137e-2
        const K35: F64x2 = F64x2::from_bits(0xBF89A2F265CB6A84, 0xBC12A235126683DC); // -1.251782773163268744627397696142e-2
        const K37: F64x2 = F64x2::from_bits(0x3F7C65EC9711289A, 0xBC147B932323683D); // 6.933139972782020946059224346584e-3
        const K39: F64x2 = F64x2::from_bits(0xBF6980ABBC25DFED, 0xBC0374FF51B89672); // -3.113112850027014096722356271323e-3
        const K41: F64x2 = F64x2::from_bits(0x3F51ABA0D832F840, 0xBBDE640F1D3AD41C); // 1.078517031416051349183701491189e-3
        const K43: F64x2 = F64x2::from_bits(0xBF319D51BE7DBE7A, 0x3BD4A9E83110A63C); // -2.687763758250982140852720885938e-4
        const K45: F64x2 = F64x2::from_bits(0x3F066187B926004B, 0x3B57E8B70AE73E4B); // 4.268832630582914715025741411327e-5
        const K47: F64x2 = F64x2::from_bits(0xBECB236D2148CA77, 0x3B54350CDAE4FFCD); // -3.235147469271306583022107945664e-6

        x + horner!(
            x3,
            x2,
            [
                K3, K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31, K33, K35,
                K37, K39, K41, K43, K45, K47
            ]
        )
    }

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

    let mut edelta = 0;
    let nexp = n.exponent();
    let dexp = d.exponent();
    if off == 0.0 && nexp - dexp < -960 {
        // Avoid subnormals in `F64x2`, which break accuracy.
        n *= f64::exp2i_fast(106);
        edelta = -106;
    }

    let z = F64x2::div11(n, d);
    if off == 0.0 && z.hi().exponent() <= -32 {
        // atan2(y, x) ~= y/x = n/d
        (z, edelta)
    } else {
        // atan2(y, x) = atan(n/d) + off * π/2
        let z2 = z * z;
        let z3 = z2 * z;
        let r = atan_poly(z, z2, z3) + off * FRAC_PI_2;
        (r, 0)
    }
}
