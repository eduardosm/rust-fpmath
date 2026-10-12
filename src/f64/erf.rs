//! Error functions for `f64`, and the parts that `f32` shares.
//!
//! `erf(-x) = -erf(x)` and `erfc(-x) = 2 - erfc(x)`, so they are calculated
//! from `|x|`.
//!
//! For `|x| < POLY_MAX`, `erf(x) = x * P(x^2)`, where `P` is approximated
//! with a polynomial (see `erf_poly`), and `erfc(x) = 1 - erf(x)`.
//!
//! Otherwise,
//!
//! `erfc(|x|) = exp(-x^2) * erfcx(|x|)`
//!
//! and `erf(|x|) = 1 - erfc(|x|)`, where `erfcx` is the scaled complementary
//! error function, which decreases slowly (it is about `1 / (sqrt(π) * x)`
//! when `x` is large). Each binade is split in four intervals of the same
//! length, and `erfcx` is approximated with a polynomial in each of them (see
//! `ERFCX_TBL`). `x^2` is calculated as a sum of two `f64`, and its
//! exponential is evaluated like in `exp` (see `super::exp`).
//!
//! The polynomials, the exponential and their products are calculated as
//! sums of two `f64`, so the relative error before the final rounding is
//! about 2^-67. When `|x| >= PLAIN_MIN`, `erfc(|x|)` is so small that plain
//! `f64` arithmetic is accurate enough to subtract it from one or two.
//!
//! `erfcx` uses the same approximations:
//!
//! * `erfcx(x) = exp(x^2) * (1 - erf(x))` for `|x| < POLY_MAX`.
//! * The polynomials of `ERFCX_TBL` for `POLY_MAX <= x < ASYMPT_MIN`.
//! * `erfcx(x) = Q(1 / x^2) / x` for `x >= ASYMPT_MIN`, where `Q` is
//!   approximated with a polynomial (see `erfcx_asympt`), as in the
//!   asymptotic series of `erfcx`.
//! * `erfcx(x) = 2 * exp(x^2) - erfcx(|x|)` for `x <= -POLY_MAX`, where the
//!   second term is small (plain `f64` arithmetic is accurate enough for it
//!   when `|x| >= PLAIN_MIN`) and negligible when `|x| >= ERFCX_EXP_MIN`.

use super::exp::{Reduced, exp_poly};
use super::f64x2::{F64x2, Split};
use super::{mul_tiny, recip_parts, square_parts};
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_2_SQRT_PI
const FRAC_2_SQRT_PI: F64x2 = F64x2::from_bits(0x3FF20DD750429B6D, 0x3C71AE3A914FED80); // 1.128379167095512573896158903122e0

/// `erf(x)` is calculated as `x * P(x^2)` for `|x| < POLY_MAX`, and from
/// `erfcx(|x|)` otherwise.
const POLY_MAX: f64 = 0.5;

/// `erfc(x) < 2^-55` for `x >= SATURATION_MIN`, so `erf(x)` rounds to one and
/// `erfc(-x)` rounds to two.
const SATURATION_MIN: f64 = 6.0;

/// `erfc(x) < 2^-20` for `x >= PLAIN_MIN`, so `1 - erfc(x)` and
/// `2 - erfc(x)` do not need `erfc(x)` as a sum of two `f64`: the error of
/// `erfc_plain` is less than 2^-70 relative to them.
///
/// Likewise, `erfcx(x) < 2^-20 * exp(x^2)`, so the error of `erfcx_plain` is
/// less than 2^-72 relative to `2 * exp(x^2) - erfcx(x)`.
const PLAIN_MIN: f64 = 3.5;

/// `erfc(x) < 2^-1076` for `x >= ERFC_ZERO_MIN`, so it rounds to zero.
const ERFC_ZERO_MIN: f64 = 27.25;

/// `erfcx(x)` is calculated with `ERFCX_TBL` for
/// `POLY_MAX <= x < ASYMPT_MIN`, and as `Q(1 / x^2) / x` otherwise.
const ASYMPT_MIN: f64 = 28.0;

/// `erfcx(x) < 2^-74 * exp(x^2)` for `x >= ERFCX_EXP_MIN`, so it is
/// negligible in `erfcx(-x) = 2 * exp(x^2) - erfcx(x)`.
const ERFCX_EXP_MIN: f64 = 7.0;

/// `2 * exp(x^2) > 2^1024` for `x > ~26.63`, so `erfcx(-x)` overflows for
/// `x >= ERFCX_INF_MIN`.
const ERFCX_INF_MIN: f64 = 26.75;

/// Approximation of `erfcx` in one of the intervals of `ERFCX_TBL`:
///
/// `erfcx(c + d) ~= K0 + K1 * d + ... + K16 * d^16`
///
/// where `c` is the center of the interval.
struct ErfcxPoly {
    /// `K0, ..., K4`
    k_hi: [F64x2; 5],
    /// `K5, ..., K16`
    k_lo: [f64; 12],
}

/// The intervals of `ERFCX_TBL` are identified by the exponent and the two
/// leading bits of the mantissa of their values: the index of the interval of
/// `x` is `(x.to_bits() >> ERFCX_TBL_SHIFT) - ERFCX_TBL_FIRST`.
const ERFCX_TBL_SHIFT: u32 = 50;

/// `x.to_bits() >> ERFCX_TBL_SHIFT` for the first interval of `ERFCX_TBL`,
/// which starts at 0.5.
const ERFCX_TBL_FIRST: u64 = 0.5f64.to_bits() >> ERFCX_TBL_SHIFT;

// GENERATE: erfcx_table 2 0.5 28 5 12
// ERFCX_TBL[i] has the coefficients `K0, ..., K16` such that
// erfcx(c + d) ~= K0 + K1 * d + ... + K16 * d^16, for |d| <= (b - a) / 2,
// where [a, b) is the i-th interval from 0.5 (each binade is split in
// 4 intervals of the same length) and c = (a + b) / 2
static ERFCX_TBL: [ErfcxPoly; 23] = [
    // [0.5, 0.625)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FE2B84F076E14FB, 0x3C8B8BD1ACCF8123), // 5.849986214749657776233867805440e-1
            F64x2::from_bits(0xBFDE18AB7052BE7F, 0x3C7C4C5D3382FB4B), // -4.702557179361760740698487750047e-1
            F64x2::from_bits(0x3FD482BD9FAD9ED0, 0xBC74FD6839673E0E), // 3.204797801358667359590968446183e-1
            F64x2::from_bits(0xBFC8BED65E01718D, 0x3C14F86D94007138), // -1.933238944065006900619045937073e-1
            F64x2::from_bits(0x3FBB1A22AA7A6DC0, 0x3C231C6BA9E3D5D5), // 1.058675447661050488996376899460e-1
        ],
        k_lo: [
            f64::from_bits(0xBFAB65943018044B), // -5.350936019022664e-2
            f64::from_bits(0x3F99DCC15144E60F), // 2.525617655303419e-2
            f64::from_bits(0xBF86FF6B38454D3F), // -1.122936025118426e-2
            f64::from_bits(0x3F7364EB29716855), // 4.734915352935761e-3
            f64::from_bits(0xBF5F301024A608B1), // -1.9035489700347903e-3
            f64::from_bits(0x3F480374A0767E3C), // 7.328338114582287e-4
            f64::from_bits(0xBF31C5283FCEE218), // -2.71150900342161e-4
            f64::from_bits(0x3F195AAC68C9A4CE), // 9.67185716917063e-5
            f64::from_bits(0xBF017B925528FDFD), // -3.334560701206619e-5
            f64::from_bits(0x3EE75B54AC0A1E55), // 1.1137371417348883e-5
            f64::from_bits(0xBECE5690F0E2C6E9), // -3.6165892033029055e-6
            f64::from_bits(0x3EB31EE64DFD15CF), // 1.1396826308946616e-6
        ],
    },
    // [0.625, 0.75)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FE0FCE4E96DD619, 0xBC67210607953137), // 5.308708724175544746040456631292e-1
            F64x2::from_bits(0xBFD97FE7BF1C60F0, 0x3C6CD62292F307BB), // -3.984317175213751713155961163164e-1
            F64x2::from_bits(0x3FD071DA7F78298D, 0xBC3A6E721470610F), // 2.569490666216090443245733331679e-1
            F64x2::from_bits(0xBFC2ECD75F4D05D4, 0x3C58C0DBAF309712), // -1.478528228126793022283013308704e-1
            F64x2::from_bits(0x3FB3E0E0ED6B5F18, 0x3C55347D3F941707), // 7.765012546894601202130805710858e-2
        ],
        k_lo: [
            f64::from_bits(0xBFA358DCE2B361D2), // -3.7787344621111565e-2
            f64::from_bits(0x3F91A31BFF52474B), // 1.7223775347310604e-2
            f64::from_bits(0xBF7E5D3E4FD95E6F), // -7.413142591381579e-3
            f64::from_bits(0x3F68D62A9331D61F), // 3.0318099539339417e-3
            f64::from_bits(0xBF5366C6D2E3C18B), // -1.1841718329002905e-3
            f64::from_bits(0x3F3D1156C45277FF), // 4.435383637629852e-4
            f64::from_bits(0xBF24F41153490674), // -1.5986167422181766e-4
            f64::from_bits(0x3F0D2741150C8C5E), // 5.560557713681838e-5
            f64::from_bits(0xBEF39F289DCB2FF8), // -1.871272421531675e-5
            f64::from_bits(0x3ED99C0B346876CF), // 6.105791561709323e-6
            f64::from_bits(0xBEC042802FCB9858), // -1.9383154480244495e-6
            f64::from_bits(0x3EA40D84E1EF22FE), // 5.976202749276064e-7
        ],
    },
    // [0.75, 0.875)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FDF0723FF5ACDF9, 0xBC57B669FFA31576), // 4.848108285616202166242079342632e-1
            F64x2::from_bits(0xBFD5CBC2C216DF00, 0x3C695F54BDABA6E5), // -3.405615706828797218818210099426e-1
            F64x2::from_bits(0x3FCAA32B83507192, 0x3C51847FB493E741), // 2.081045523817804425952283636875e-1
            F64x2::from_bits(0xBFBD43ED8CA5D7A1, 0xBC440793AFC061DB), // -1.143177479151220748487986591122e-1
            F64x2::from_bits(0x3FAD7F26045A23F1, 0x3C45616BB37E3D79), // 5.761069110037187839028971568660e-2
        ],
        k_lo: [
            f64::from_bits(0xBF9BA6D691685B3F), // -2.700362455842797e-2
            f64::from_bits(0x3F8859D3CC5FA92B), // 1.1890082048883052e-2
            f64::from_bits(0xBF744BD2E2B89098), // -4.955123683917283e-3
            f64::from_bits(0x3F601B062044AE6D), // 1.966011013925065e-3
            f64::from_bits(0xBF487347D1168039), // -7.461643855784023e-4
            f64::from_bits(0x3F31D29276135FD2), // 2.7195049012851164e-4
            f64::from_bits(0xBF190857BBFACB1C), // -9.549174773830184e-5
            f64::from_bits(0x3F00FBD62FEEEE33), // 3.2393907523933437e-5
            f64::from_bits(0xBEE651445A0C50B5), // -1.064178906199534e-5
            f64::from_bits(0x3ECC7551472D451E), // 3.3924903219204047e-6
            f64::from_bits(0xBEB1AA555348A270), // -1.052937773535691e-6
            f64::from_bits(0x3E95512F0779BB68), // 3.176498988071293e-7
        ],
    },
    // [0.875, 1)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FDC7F81382721EF, 0xBC13F947EE724DFD), // 4.452822731368816786059247755797e-1
            F64x2::from_bits(0xBFD2C84AF7C10E14, 0xBC64881059988BD0), // -2.934749049638594265100499489089e-1
            F64x2::from_bits(0x3FC5C775DFC44978, 0x3C58414CD91ADAB0), // 1.701495497332634662527529484787e-1
            F64x2::from_bits(0xBFB6DCC9E7A2C9B5, 0x3C423C1F1F241619), // -8.930646805928328459872938165619e-2
            F64x2::from_bits(0x3FA61FEE765FF5D7, 0xBC36650921E5DCFC), // 4.321236796384269347072207237989e-2
        ],
        k_lo: [
            f64::from_bits(0xBF93FC8366A2E3F3), // -1.9517949237272302e-2
            f64::from_bits(0x3F8102167DAF7951), // 8.304763517966636e-3
            f64::from_bits(0xBF6B75CB1ACA9956), // -3.3520666969081665e-3
            f64::from_bits(0x3F5524F5C6CFFAC2), // 1.2905502474038076e-3
            f64::from_bits(0xBF3F329C42D74AD5), // -4.760390755482018e-4
            f64::from_bits(0x3F2621C1D83C073C), // 1.6885272281546594e-4
            f64::from_bits(0xBF0E49DF1F6C78B2), // -5.777084509460527e-5
            f64::from_bits(0x3EF40B4211408EF0), // 1.9115425928486107e-5
            f64::from_bits(0xBED9B6E11CCDD0BD), // -6.130784140679656e-6
            f64::from_bits(0x3EC005051986A471), // 1.9096862155688163e-6
            f64::from_bits(0xBEA3724288C42CC9), // -5.795457283254073e-7
            f64::from_bits(0x3E86F670A4E1AB8F), // 1.7108512157261125e-7
        ],
    },
    // [1, 1.25)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FD9531E09B149B5, 0xBC5AA513235EB798), // 3.956980795529959014044473809124e-1
            F64x2::from_bits(0xBFCE78B356770FBB, 0x3C6EA9D5564683DB), // -2.380584881012717957361522782036e-1
            F64x2::from_bits(0x3FC05E72521CA1B8, 0xBC369BCC7F07CD04), // 1.278822804390651312012760872249e-1
            F64x2::from_bits(0xBFB01343A2C92265, 0x3C4E9FAB05861DB1), // -6.279394840488234875653269669512e-2
            F64x2::from_bits(0x3F9D4E711A2DF97D, 0x3C1C96A68C8CA543), // 2.861954424178624442507474173255e-2
        ],
        k_lo: [
            f64::from_bits(0xBF8910A5D7C0A71F), // -1.223878445314913e-2
            f64::from_bits(0x3F7446C5166CCF50), // 4.950303910664491e-3
            f64::from_bits(0xBF5F38C6D0510616), // -1.9056264439004703e-3
            f64::from_bits(0x3F46FD9A57AC0B52), // 7.016185403191189e-4
            f64::from_bits(0xBF3041E38D545D68), // -2.48067908004987e-4
            f64::from_bits(0x3F162743C046275C), // 8.450842876164599e-5
            f64::from_bits(0xBEFD2B3006D37AAF), // -2.781744153856441e-5
            f64::from_bits(0x3EE2997DADF71F37), // 8.868967994582137e-6
            f64::from_bits(0xBEC705E4672CCA62), // -2.7445574287907215e-6
            f64::from_bits(0x3EABB65F6AA8A15D), // 8.258936835025427e-7
            f64::from_bits(0xBE9054ABAB10D03C), // -2.433470554536376e-7
            f64::from_bits(0x3E72A7EF2C62E244), // 6.949899064829014e-8
        ],
    },
    // [1.25, 1.5)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FD5F88F52F3C76B, 0xBC6B7EB97A02D456), // 3.432958898621253798874683250608e-1
            F64x2::from_bits(0xBFC797A639D8129D, 0xBC5DF1E663E51527), // -1.843154699746677792056210038432e-1
            F64x2::from_bits(0x3FB701342CBCEA7B, 0xBC502A0FBD328502), // 8.986211864695718347973944710530e-2
            F64x2::from_bits(0xBFA4BCDB9D9083C2, 0xBC403A5A004BB7E2), // -4.050337122340110128066931895253e-2
            F64x2::from_bits(0x3F917EBA60D31FCA, 0x3C395BA5A0D4BBD2), // 1.708499160739033460941027644323e-2
        ],
        k_lo: [
            f64::from_bits(0xBF7BDF24BCCAC617), // -6.804603105295757e-3
            f64::from_bits(0x3F651AB9FFCE7487), // 2.57622077920289e-3
            f64::from_bits(0xBF4E8AE68B41E94E), // -9.32085581111944e-4
            f64::from_bits(0x3F3535F57FDF98C7), // 3.236507762934934e-4
            f64::from_bits(0xBF1C5FA6B09B5A20), // -1.0823683637837438e-4
            f64::from_bits(0x3F0254ED1EA8F86B), // 3.49650252542279e-5
            f64::from_bits(0xBEE6F06273424543), // -1.0938168635984464e-5
            f64::from_bits(0x3ECBDB736FE93F52), // 3.3208406187492474e-6
            f64::from_bits(0xBEB0725A3CA44BF8), // -9.80299079748847e-7
            f64::from_bits(0x3E92EA0695D99E53), // 2.8184299147866206e-7
            f64::from_bits(0xBE7551DBB074B2E0), // -7.942228929099891e-8
            f64::from_bits(0x3E575670D36C8DEE), // 2.173488874568589e-8
        ],
    },
    // [1.5, 1.75)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FD3583F6644327B, 0xBC788EB8EBFDCD23), // 3.022612093634859174056627734060e-1
            F64x2::from_bits(0xBFC2B11E6959934C, 0x3C6D03D8DF7F2509), // -1.460302366641833423277548878649e-1
            F64x2::from_bits(0x3FB0A15AC2ADAB35, 0xBC2F7316CCF36753), // 6.496207478418798612306108040128e-2
            F64x2::from_bits(0xBF9BA018E6428103, 0x3C236C4EBE38CB98), // -2.697791009325190991852560596837e-2
            F64x2::from_bits(0x3F85A142948A9B2F, 0x3C101CF0B44BB169), // 1.056148544132681625273004611834e-2
        ],
        k_lo: [
            f64::from_bits(0xBF7014EAE28304AE), // -3.926198500438334e-3
            f64::from_bits(0x3F56D609F6AB13B0), // 1.3938042927048415e-3
            f64::from_bits(0xBF3F1B43D3AB833F), // -4.746475785122779e-4
            f64::from_bits(0x3F2465ECD15ACCD8), // 1.5562549440559822e-4
            f64::from_bits(0xBF09D62282C94D44), // -4.92791444669746e-5
            f64::from_bits(0x3EEFAFC8F3E83950), // 1.510937692919255e-5
            f64::from_bits(0xBED2DB3B7752CE65), // -4.495710404106913e-6
            f64::from_bits(0x3EB5D23634396CFD), // 1.300641273702465e-6
            f64::from_bits(0xBE98982626F1E583), // -3.6648410790088477e-7
            f64::from_bits(0x3E7B09F6AA0B0CB9), // 1.0072782653373366e-7
            f64::from_bits(0xBE5D2A35BB54D0DC), // -2.7161913348023422e-8
            f64::from_bits(0x3E3E9ACF358B445B), // 7.125717648493065e-9
        ],
    },
    // [1.75, 2)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FD13E5743B60480, 0x3C7CA1DFCA5D5311), // 2.694299851646704517460504372285e-1
            F64x2::from_bits(0xBFBE36580C7F734A, 0xBC593CCD69BCD773), // -1.180167227279983798484697629592e-1
            F64x2::from_bits(0x3FA8A6EFEED233AD, 0xBC499503F9A79F5A), // 4.814863004967348953016963132498e-2
            F64x2::from_bits(0xBF92EF92F6F10797, 0x3C1DA540EBFA5EAA), // -1.849202758990705798626950963525e-2
            F64x2::from_bits(0x3F7B99589D40F23D, 0xBC1A891EA40690B8), // 6.738039159298877902958059441302e-3
        ],
        k_lo: [
            f64::from_bits(0xBF633237C3EEACEE), // -2.3432816664886646e-3
            f64::from_bits(0x3F499B60E42DD5A3), // 7.814620115442105e-4
            f64::from_bits(0xBF3070E0CB5E266B), // -2.508686842409349e-4
            f64::from_bits(0x3F14631C4B0B2353), // 7.777080714811467e-5
            f64::from_bits(0xBEF87A61E43B99E5), // -2.334409351947312e-5
            f64::from_bits(0x3EDC8594802F87ED), // 6.80012635975896e-6
            f64::from_bits(0xBEC0286353E73DC6), // -1.9261557604244667e-6
            f64::from_bits(0x3EA1D4F4862A6D33), // 5.314307254979567e-7
            f64::from_bits(0xBE833295A4A14AC0), // -1.4303323531363736e-7
            f64::from_bits(0x3E64307CFAEF5765), // 3.760570110144271e-8
            f64::from_bits(0xBE44D9FC4AD94ECB), // -9.709739091891504e-9
            f64::from_bits(0x3E24FC834E01C28C), // 2.4431360129894747e-9
        ],
    },
    // [2, 2.5)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FCD94446D627932, 0xBC5A8198A876279A), // 2.310872587303918699574998886172e-1
            F64x2::from_bits(0xBFB6A70D2BB37411, 0x3C1FFE90ACFEC198), // -8.848650280874915908740188038690e-2
            F64x2::from_bits(0x3FA0615670E25A7B, 0xBC38073D20EA5B75), // 3.199262741070626201083583375206e-2
            F64x2::from_bits(0xBF86883F9919A17A, 0x3C27B3315BB86BB1), // -1.100206075644004638112983234591e-2
            F64x2::from_bits(0x3F6DA595561F7D33, 0x3BF6816AAA51ECB4), // 3.618995354358078831498655937173e-3
        ],
        k_lo: [
            f64::from_bits(0xBF52BD251BB2F023), // -1.1437284836537463e-3
            f64::from_bits(0x3F36D7743D3B280F), // 3.485354220457157e-4
            f64::from_bits(0xBF1AED7EBC9A0859), // -1.0272108115752641e-4
            f64::from_bits(0x3EFEC773CC926889), // 2.9353247360399313e-5
            f64::from_bits(0xBEE117A66607E05B), // -8.150283236738208e-6
            f64::from_bits(0x3EC27AF428B6491F), // 2.2030220116566087e-6
            f64::from_bits(0xBEA37B9AD1161995), // -5.806336124548889e-7
            f64::from_bits(0x3E840E792D62C52C), // 1.4943284182758858e-7
            f64::from_bits(0xBE642F69517E939B), // -3.7597866319676296e-8
            f64::from_bits(0x3E43E31D62FA2EF5), // 9.260683858016677e-9
            f64::from_bits(0xBE237B72E46DBE35), // -2.268029128592888e-9
            f64::from_bits(0x3E026A3A66E80E4A), // 5.359456887081352e-10
        ],
    },
    // [2.5, 3)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FC8C9EB68FF27D7, 0xBC6BB4E763CAA5C2), // 1.936620962790686786026834799464e-1
            F64x2::from_bits(0xBFB0305781330099, 0x3C4102498310AF5C), // -6.323763756063484158139872312091e-2
            F64x2::from_bits(0x3F943B98BAC83823, 0xBBBAA96D8F4F8FFE), // 1.975859298732286425383481973838e-2
            F64x2::from_bits(0xBF784E9AB30E6AB3, 0xBC10C66B97BE4BD6), // -5.934337896997976589701280717796e-3
            F64x2::from_bits(0x3F5C2C72FD72763E, 0x3BC98065CFEAB2DA), // 1.719581885289214317039741710901e-3
        ],
        k_lo: [
            f64::from_bits(0xBF3F99E41ECB0901), // -4.8219508498105473e-4
            f64::from_bits(0x3F2131BB16125574), // 1.311818005304378e-4
            f64::from_bits(0xBF02312B2596801D), // -3.4698609577832525e-5
            f64::from_bits(0x3EE2BFB5B0EB951D), // 8.940156047863298e-6
            f64::from_bits(0xBEC2DA329C266419), // -2.2473734315295592e-6
            f64::from_bits(0x3EA2856FAB0D0952), // 5.519758215849241e-7
            f64::from_bits(0xBE81CCF9F777F905), // -1.3262546948162578e-7
            f64::from_bits(0x3E60C1602337A615), // 3.120931460125242e-8
            f64::from_bits(0xBE3EEBE86D524963), // -7.199476311245849e-9
            f64::from_bits(0x3E1C0007A4425CD3), // 1.6298212927147063e-9
            f64::from_bits(0xBDF936CC8B2B29DF), // -3.6691284313658504e-10
            f64::from_bits(0x3DD602E8A4F4F2E0), // 8.007686982184464e-11
        ],
    },
    // [3, 3.5)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FC54A7A08D4BB45, 0xBC26A0D9133F2E3E), // 1.663353484268218767633828671151e-1
            F64x2::from_bits(0xBFA82A8522B868A1, 0x3C3B91959FCDC7F0), // -4.719940232117037493417010094783e-2
            F64x2::from_bits(0x3F8A7EDDC9EE6425, 0x3C2C91725DD4CCF0), // 1.293729088301815822732957754730e-2
            F64x2::from_bits(0xBF6C24B49C47A2C4, 0x3BF535AE1281C66D), // -3.435471300907573797026784789932e-3
            F64x2::from_bits(0x3F4D085857A17F33, 0xBBE6688138EEBB42), // 8.860045775342716936912224751621e-4
        ],
        k_lo: [
            f64::from_bits(0xBF2D25EBBA1C4910), // -2.2238256956847629e-4
            f64::from_bits(0x3F0C882F0238146E), // 5.442040881224122e-5
            f64::from_bits(0xBEEB45D025FA2D50), // -1.3004640265343537e-5
            f64::from_bits(0x3EC97DD78D735668), // 3.038831987471278e-6
            f64::from_bits(0xBEA753CAB56B91C1), // -6.952080678605071e-7
            f64::from_bits(0x3E84EC091FE16C61), // 1.558811532656179e-7
            f64::from_bits(0xBE6268C3EE2E1B28), // -3.428988095467454e-8
            f64::from_bits(0x3E3FCF8B34D9B94D), // 7.4065099171861065e-9
            f64::from_bits(0xBE1B01DE0B74878A), // -1.5720314336914217e-9
            f64::from_bits(0x3DF68CE1A929262B), // 3.2815032042424427e-10
            f64::from_bits(0xBDD2BAAC006AF1A0), // -6.813638551327514e-11
            f64::from_bits(0x3DAE47E62F09A7C4), // 1.3770139082749091e-11
        ],
    },
    // [3.5, 4)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FC2A2AF19C14930, 0xBC6FA04A06A35214), // 1.455897212750385390456688240194e-1
            F64x2::from_bits(0xBFA2AA6503ACDA11, 0xBC11D40A67D817A5), // -3.645625753272353105364269288468e-2
            F64x2::from_bits(0x3F822F0664F3CBF9, 0xBC2B72342C905889), // 8.878755527325297594508622872043e-3
            F64x2::from_bits(0xBF61434AE05873AB, 0x3C05C3B37FEFE33C), // -2.107282870169110049513352806108e-3
            F64x2::from_bits(0x3F3FFF032A0DF889, 0x3BCEC21A3DB4B438), // 4.882223820955674544591085994662e-4
        ],
        k_lo: [
            f64::from_bits(0xBF1CFCDEA1B1F551), // -1.1057957492429284e-4
            f64::from_bits(0x3EF9B50D0D260D9C), // 2.451632537648977e-5
            f64::from_bits(0xBED65778AAD3973A), // -5.326672789273719e-6
            f64::from_bits(0x3EB30C2FB3FEC951), // 1.1353256041788716e-6
            f64::from_bits(0xBE8FE3E32B2E0D90), // -2.3760039410620524e-7
            f64::from_bits(0x3E6A3BEE3168E3C6), // 4.8864825231667635e-8
            f64::from_bits(0xBE4539512C6BDD95), // -9.883146202020738e-9
            f64::from_bits(0x3E20E5DB47DD1A77), // 1.9671718133862722e-9
            f64::from_bits(0xBDFA7EE429D01945), // -3.8556272915081997e-10
            f64::from_bits(0x3DD477DD6335E080), // 7.446295733020398e-11
            f64::from_bits(0xBDAF7771A18270A2), // -1.4309342814104772e-11
            f64::from_bits(0x3D879F181E4415EC), // 2.685449290416871e-12
        ],
    },
    // [4, 5)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FBF5B2A049CF4C6, 0xBC5FC4BBCA0E3564), // 1.224848042738414175492253205526e-1
            F64x2::from_bits(0xBF9AA3EB6A946F7E, 0xBC328D280F052C42), // -2.601592863093981595293792325973e-2
            F64x2::from_bits(0x3F762C12CB5F7577, 0xBC0490541E115AF7), // 5.413125434612245760178235079686e-3
            F64x2::from_bits(0xBF5218ED930B23B1, 0xBBE86C9D4EA5ED04), // -1.104576116789806725183560884503e-3
            F64x2::from_bits(0x3F2D00785F2ED4FB, 0x3BA0B039C543823B), // 2.212664545290578316474489980848e-4
        ],
        k_lo: [
            f64::from_bits(0xBF06D54B133D63CD), // -4.3550828563616516e-5
            f64::from_bits(0x3EE1AD6A32DEBC43), // 8.429241997591252e-6
            f64::from_bits(0xBEBAEF8EA4AFC42C), // -1.6054970213253735e-6
            f64::from_bits(0x3E94354FB1F1665C), // 3.0112635047083114e-7
            f64::from_bits(0xBE6DE095362BD0D0), // -5.5650764682186324e-8
            f64::from_bits(0x3E45C64B2FBEEFE4), // 1.0139581167629917e-8
            f64::from_bits(0xBE1F4E95387C3F86), // -1.8223051695595026e-9
            f64::from_bits(0x3DF635E663D29BD9), // 3.232059930652087e-10
            f64::from_bits(0xBDCF1A12733E9338), // -5.6573924716406816e-11
            f64::from_bits(0x3DA584AF1E63D813), // 9.785388601572792e-12
            f64::from_bits(0xBD7E53C2592DF126), // -1.7239008383619541e-12
            f64::from_bits(0x3D5475539C8BA1CD), // 2.907300293533039e-13
        ],
    },
    // [5, 6)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FB9D8A8F2284F2C, 0xBC55B027800B2C88), // 1.009622183994990882327984899818e-1
            F64x2::from_bits(0xBF9238CA71B93FC3, 0xBC352A64DB942BD6), // -1.779476470102260333536374106409e-2
            F64x2::from_bits(0x3F695252B932EFE3, 0xBC0950AAF0EB4910), // 3.091012543874769888234511197663e-3
            F64x2::from_bits(0xBF415976DDDA3CA1, 0x3BCF6092E1CA9030), // -5.294638064742459689768576273675e-4
            F64x2::from_bits(0x3F1774F4826DC857, 0xBBA595D21E59EDD7), // 8.948080413320853571467205761815e-5
        ],
        k_lo: [
            f64::from_bits(0xBEEF4E46D179E36D), // -1.4927753496639482e-5
            f64::from_bits(0x3EC4A17E19BFD2EE), // 2.4593866338968857e-6
            f64::from_bits(0xBE9ADD7AC9CB40ED), // -4.0032200291935943e-7
            f64::from_bits(0x3E7149CD78261DFA), // 6.440390446480369e-8
            f64::from_bits(0xBE460001E6372DA1), // -1.024456181596245e-8
            f64::from_bits(0x3E1BB09BA2BA3F08), // 1.6117628437711955e-9
            f64::from_bits(0xBDF13D9DF45D7377), // -2.508850803886075e-10
            f64::from_bits(0x3DC53F6E7399796E), // 3.864948504463346e-11
            f64::from_bits(0xBD99EAD3BAC4C2B8), // -5.8929102207839645e-12
            f64::from_bits(0x3D6F4FDF7FBA8EB5), // 8.899406814517509e-13
            f64::from_bits(0xBD432D5B23B3B571), // -1.3626200205727173e-13
            f64::from_bits(0x3D16B65CE03CD69E), // 2.0172623282199008e-14
        ],
    },
    // [6, 7)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FB5F75C42E97171, 0xBC5A1EAAA81F223D), // 8.580567010489460177788758580559e-2
            F64x2::from_bits(0xBF8A6E2CF277A0CB, 0xBC223CA411E0C01B), // -1.290545573188275078361929767995e-2
            F64x2::from_bits(0x3F5F75EF7A0FB2DB, 0x3BD00E1BAB79956A), // 1.920207847656721684355770186069e-3
            F64x2::from_bits(0xBF328787F0E7F4C6, 0xBBD509FE46B4C49C), // -2.827364814093732237279663291781e-4
            F64x2::from_bits(0x3F059B29296C518D, 0xBBACE0126D6B2000), // 4.121035924789786568871643552561e-5
        ],
        k_lo: [
            f64::from_bits(0xBED8F23FECCB7390), // -5.947658519214828e-6
            f64::from_bits(0x3EAC871A15EAF7B8), // 8.501929576671372e-7
            f64::from_bits(0xBE8028F46472CEC5), // -1.2040122696553526e-7
            f64::from_bits(0x3E52246795855E33), // 1.6896245598250682e-8
            f64::from_bits(0xBE2430034610DAFD), // -2.350140124631205e-9
            f64::from_bits(0x3DF6450BCE96DF90), // 3.240669526038173e-10
            f64::from_bits(0xBDC85C14872DF4F3), // -4.4310014769704965e-11
            f64::from_bits(0x3D9A6D2AB053A821), // 6.00867511571686e-12
            f64::from_bits(0xBD6C6F05E5A22C30), // -8.081338970432546e-13
            f64::from_bits(0x3D3E5B6CDCA053C9), // 1.0785019046568929e-13
            f64::from_bits(0xBD1060DF5051D315), // -1.4546948088231008e-14
            f64::from_bits(0x3CE135BB16284E80), // 1.910681164867422e-15
        ],
    },
    // [7, 8)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FB31742F4D8D4D3, 0xBC4F82BC955F2586), // 7.457369306287668300512924608336e-2
            F64x2::from_bits(0xBF8404455BA9F3BE, 0xBC286B84E93C1858), // -9.773771152362328819220103535939e-3
            F64x2::from_bits(0x3F54D07BBA601411, 0x3BF67536FD4EF2D0), // 1.270409420159216860977659202785e-3
            F64x2::from_bits(0xBF25783A27FAD25A, 0xBBBF4EEF4ECFFAAB), // -1.638003341121349079459088215095e-4
            F64x2::from_bits(0x3EF5F8A69C9DDB98, 0xBB75E9A32998741F), // 2.095345715910252577082007448707e-5
        ],
        k_lo: [
            f64::from_bits(0xBEC64FCBB7153F92), // -2.6597621675463845e-6
            f64::from_bits(0x3E967CA2BF7EAA56), // 3.3508030083487716e-7
            f64::from_bits(0xBE667F14D3525844), // -4.1902831795688136e-8
            f64::from_bits(0x3E3657F396188984), // 5.202265591861659e-9
            f64::from_bits(0xBE0608EA24D96E43), // -6.412977455484323e-10
            f64::from_bits(0x3DD594675DC0CC99), // 7.850649942482038e-11
            f64::from_bits(0xBDA4FD83324F095E), // -9.545275656819998e-12
            f64::from_bits(0x3D7447DDE3A6C94C), // 1.1528260024232476e-12
            f64::from_bits(0xBD43770C54E0CE56), // -1.3830737028257997e-13
            f64::from_bits(0x3D1290539BA9E5BF), // 1.6487945016203173e-14
            f64::from_bits(0xBCE1DA7BC0D0875A), // -1.9821312173482885e-15
            f64::from_bits(0x3CB0CEB47D897884), // 2.33250107626027e-16
        ],
    },
    // [8, 10)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FAFE6CC4A0F9538, 0x3C4C95BE413B3EE5), // 6.230772403777468414653727515049e-2
            F64x2::from_bits(0xBF7C046699D77D51, 0xBC174C1D9F527928), // -6.840134415568259258090487435095e-3
            F64x2::from_bits(0x3F4876373F4A0F3A, 0x3BE5EE3B52822C68), // 7.465142976603508201693157684323e-4
            F64x2::from_bits(0xBF153C131B5CB2A0, 0xBBBF7707BCD5E770), // -8.100382441673460346017198209807e-5
            F64x2::from_bits(0x3EE254381038C8EC, 0x3B7CCD32A7041276), // 8.739938954869780841375921099871e-6
        ],
        k_lo: [
            f64::from_bits(0xBEAF773505C3477D), // -9.377495291623631e-7
            f64::from_bits(0x3E7ADC5F7F272D9B), // 1.000643974686887e-7
            f64::from_bits(0xBE46CE6631B7CACA), // -1.0619986271441175e-8
            f64::from_bits(0x3E1342C874A5A207), // 1.1211302603196651e-9
            f64::from_bits(0xBDE02E7B6BD5B90E), // -1.177364229698499e-10
            f64::from_bits(0x3DAB0C8E0548B603), // 1.2300480222165535e-11
            f64::from_bits(0xBD767E2FECA21E0E), // -1.2785744029047239e-12
            f64::from_bits(0x3D429C461ED1F47E), // 1.3223516470089924e-13
            f64::from_bits(0xBD0E9E6D93C931B2), // -1.3597505129296806e-14
            f64::from_bits(0x3CD91645A3B6D3D2), // 1.3926082572984909e-15
            f64::from_bits(0xBCA56294DEF0B360), // -1.483888296200864e-16
            f64::from_bits(0x3C7156C1E748545E), // 1.5039095170416835e-17
        ],
    },
    // [10, 12)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FAA273C57E476CE, 0xBC46D026E20C4A0C), // 5.108059475808844370998557847796e-2
            F64x2::from_bits(0xBF72DDD3D589BFB0, 0xBC1AA91B35670502), // -4.606082417566812276462005460898e-3
            F64x2::from_bits(0x3F3B1C8923879D89, 0x3BABD71B710DFFFE), // 4.136881648535086687481269513160e-4
            F64x2::from_bits(0xBF03672F098120BB, 0xBB8FDC789ADB7964), // -3.700840278547794750202686605296e-5
            F64x2::from_bits(0x3ECBAA1D9D8421F9, 0x3B69F68D33EA5D79), // 3.297867106625626872866699106627e-6
        ],
        k_lo: [
            f64::from_bits(0xBE93A555D7A845A1), // -2.9274584503841124e-7
            f64::from_bits(0x3E5BCBEE2A203594), // 2.5887603734332384e-8
            f64::from_bits(0xBE23972862992DFA), // -2.2806297031327374e-9
            f64::from_bits(0x3DEB82D452160A2A), // 2.0016925013619504e-10
            f64::from_bits(0xBDB33EEF565238E5), // -1.7503989052537914e-11
            f64::from_bits(0x3D7AD45237E9C3EA), // 1.5250736597287363e-12
            f64::from_bits(0xBD42A21661244A87), // -1.3239652207246595e-13
            f64::from_bits(0x3D09CA08C99AA2F7), // 1.1452703934922252e-14
            f64::from_bits(0xBCD1C6C44E2BA197), // -9.867902540637422e-16
            f64::from_bits(0x3C986E5C122D2FB3), // 8.4762379032315e-17
            f64::from_bits(0xBC61430CDB5B46CA), // -7.486162267406202e-18
            f64::from_bits(0x3C278CCFF6BE9C02), // 6.383253375868639e-19
        ],
    },
    // [12, 14)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FA627BCC27CB4AE, 0xBC3D520E73F94CFC), // 4.327192186460969266318759881543e-2
            F64x2::from_bits(0xBF6B1BE47A913FC3, 0x3C06370120C18722), // -3.309198615660564653280505416224e-3
            F64x2::from_bits(0x3F30898F73586FCF, 0xBBCC5CCF3C4069F3), // 2.523398610223521705302759440521e-4
            F64x2::from_bits(0xBEF41E750B2E0636, 0xBB9CE2BC3412A703), // -1.918694824665762429771590984960e-5
            F64x2::from_bits(0x3EB8682D21856146, 0xBB3372525B9B53B4), // 1.454766907901527589560116233604e-6
        ],
        k_lo: [
            f64::from_bits(0xBE7D868CC70FABB1), // -1.099913775751057e-7
            f64::from_bits(0x3E41CF1FF6C645F3), // 8.292999808382074e-9
            f64::from_bits(0xBE056CB0DD5B175A), // -6.235371617574222e-10
            f64::from_bits(0x3DC9B410BB2C06FF), // 4.6754176395451367e-11
            f64::from_bits(0xBD8EC0BA0E207AB9), // -3.496193014233637e-12
            f64::from_bits(0x3D5258F4A06057A6), // 2.607334110920458e-13
            f64::from_bits(0xBD15D5827E32373C), // -1.9392507523767874e-14
            f64::from_bits(0x3CD9E9FE4A8271C7), // 1.438517994919289e-15
            f64::from_bits(0xBC9EAB26D2B8839F), // -1.0640294598157742e-16
            f64::from_bits(0x3C621AB5C9495ED9), // 7.851504315816582e-18
            f64::from_bits(0xBC25D1438A7B3328), // -5.913628105809537e-19
            f64::from_bits(0x3BE99E820F200D34), // 4.3400565605475596e-20
        ],
    },
    // [14, 16)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FA33714A024097E, 0x3C4F3B4D0AE35594), // 3.752960638850576574606117817304e-2
            F64x2::from_bits(0xBF6467F441A50CF9, 0xBBEFB34E3F16033E), // -2.490975440339601514323487676402e-3
            F64x2::from_bits(0x3F259FA2994D54EF, 0x3BBA8C22620F26D5), // 1.649747834117430312078102411198e-4
            F64x2::from_bits(0xBEE6DD369DB39094, 0x3B58E8B59DB38CCA), // -1.090245944230403080759501317799e-5
            F64x2::from_bits(0x3EA81FB2BA7EB241, 0xBB47FD9FA892A2E4), // 7.189458885912845723210627512689e-7
        ],
        k_lo: [
            f64::from_bits(0xBE696605F1EA2FB8), // -4.730844537390484e-8
            f64::from_bits(0x3E2AAF0D3BB729D6), // 3.1064026609037434e-9
            f64::from_bits(0xBDEBF9951F90513C), // -2.0354441724278e-10
            f64::from_bits(0x3DAD445B88020F1D), // 1.3309100566637327e-11
            f64::from_bits(0xBD6E8E11C9B96884), // -8.684241641707441e-13
            f64::from_bits(0x3D2FD55EDDF3CFEA), // 5.654761779816886e-14
            f64::from_bits(0xBCF08C716DC8597A), // -3.674528633971852e-15
            f64::from_bits(0x3CB12B9A8EBDECC8), // 2.382861562454864e-16
            f64::from_bits(0xBC71C6DF293CC9FF), // -1.5418953152380425e-17
            f64::from_bits(0x3C325EB7750D7C37), // 9.95838981767124e-19
            f64::from_bits(0xBBF34878CC56C2C9), // -6.533366117749624e-20
            f64::from_bits(0x3BB3D8015F3284E3), // 4.202081944503008e-21
        ],
    },
    // [16, 20)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3FA005FE090C5E5D, 0xBC44A17D59AB3112), // 3.129571781590520988566471457291e-2
            F64x2::from_bits(0xBF5C661852C5133A, 0xBBD37DF79031BD8B), // -1.733325722925018012028642483266e-3
            F64x2::from_bits(0x3F1920B4FB071805, 0xBBAE97E577A26915), // 9.585480325488566557273484122528e-5
            F64x2::from_bits(0xBED633259365787E, 0x3B718A66C10C69E4), // -5.292842891384023552870369142227e-6
            F64x2::from_bits(0x3E93955A9A783D3F, 0x3B3E0A3903EE6ACA), // 2.918156049866423355630301489108e-7
        ],
        k_lo: [
            f64::from_bits(0xBE513FDBD8E16562), // -1.606480064977618e-8
            f64::from_bits(0x3E0E578451C00686), // 8.830644301733087e-10
            f64::from_bits(0xBDCAA561A7C4CEA4), // -4.846883048655852e-11
            f64::from_bits(0x3D875D9D49719F29), // 2.656370413373736e-12
            f64::from_bits(0xBD44757D7C648A44), // -1.4536955470857706e-13
            f64::from_bits(0x3D01E3333905694C), // 7.943645895277234e-15
            f64::from_bits(0xBCBF3BAE88E4BFF8), // -4.3344677761998e-16
            f64::from_bits(0x3C7B3A53CFD4B1B1), // 2.3616387809905772e-17
            f64::from_bits(0xBC37AD1976AE7614), // -1.2834877363229665e-18
            f64::from_bits(0x3BF494A0360058D1), // 6.972968218248197e-20
            f64::from_bits(0xBBB2CCCC6B106607), // -3.981053618494178e-21
            f64::from_bits(0x3B704A432D853BA9), // 2.1559751359817595e-22
        ],
    },
    // [20, 24)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3F9A3BC12161D4EB, 0xBC2B352D341C9078), // 2.561857000587945266761757843674e-2
            F64x2::from_bits(0xBF530A253D642E9C, 0xBBF758495D32CD9D), // -1.162086836816656520979682238485e-3
            F64x2::from_bits(0x3F0B9BD9F029A792, 0x3BA6D8091DBBAF18), // 5.265959591300920593842014558469e-5
            F64x2::from_bits(0xBEC3FF35B4B9DFA9, 0x3B476ABC7E390679), // -2.383817820302660292197484442940e-6
            F64x2::from_bits(0x3E7CF017AA140A55, 0xBAFE49A9C7F4C0C6), // 1.078019331753405135007096275586e-7
        ],
        k_lo: [
            f64::from_bits(0xBE34EABFD67A0A7F), // -4.8701161780673554e-9
            f64::from_bits(0x3DEE35421D7C37DB), // 2.197924192844557e-10
            f64::from_bits(0xBDA5CA8142AA8245), // -9.909415374475114e-12
            f64::from_bits(0x3D5F6830EB5B5CFC), // 4.4632026359520693e-13
            f64::from_bits(0xBD169C477817B5E4), // -2.0082127544985923e-14
            f64::from_bits(0x3CD042EC64E30B98), // 9.026901219592225e-16
            f64::from_bits(0xBC875DFEE4FBD2B3), // -4.053557946060057e-17
            f64::from_bits(0x3C40C5B056A4A687), // 1.8184479487753113e-18
            f64::from_bits(0xBBF80A9E579067BA), // -8.145569792639596e-20
            f64::from_bits(0x3BB138E0C733828C), // 3.646938437442489e-21
            f64::from_bits(0xBB6983D683D1640A), // -1.6884405361416717e-22
            f64::from_bits(0x3B223D4065EEC630), // 7.543582737763521e-24
        ],
    },
    // [24, 28)
    ErfcxPoly {
        k_hi: [
            F64x2::from_bits(0x3F963438BF2C3EE9, 0xBC21FC7374778197), // 2.168358485056290661617299165464e-2
            F64x2::from_bits(0xBF4B49A775427EF7, 0xBBE0EC51B0B230AB), // -8.327548662414298551630410686613e-4
            F64x2::from_bits(0x3F00C15FCC6F7FEA, 0x3B7BAA2DA26BC43D), // 3.195832828573038192626162793474e-5
            F64x2::from_bits(0xBEB48FB6BC69793C, 0xBB5368BB5F85F0DB), // -1.225553874959950057044688056976e-6
            F64x2::from_bits(0x3E6936A6B39AD224, 0x3B0DA46DB5E282C7), // 4.696376838584026759147773489572e-8
        ],
        k_lo: [
            f64::from_bits(0xBE1EE543F91730D9), // -1.7983587712412274e-9
            f64::from_bits(0x3DD2EA50CFFC2CD5), // 6.881344452267718e-11
            f64::from_bits(0xBD8724F1C7A9E101), // -2.631203900482591e-12
            f64::from_bits(0x3D3C4C5DC2BF892E), // 1.0053577765944267e-13
            f64::from_bits(0xBCF14999851B17BE), // -3.838595836403411e-15
            f64::from_bits(0x3CA51B502E01F803), // 1.464570978062869e-16
            f64::from_bits(0xBC59C047F1AE896A), // -5.583879231258162e-18
            f64::from_bits(0x3C0F6514984939F4), // 2.127397481182112e-19
            f64::from_bits(0xBBC31E812C08EDBB), // -8.09727859022671e-21
            f64::from_bits(0x3B7746AB62808164), // 3.080559914787773e-22
            f64::from_bits(0xBB2D07315C0DC533), // -1.2005739494732131e-23
            f64::from_bits(0x3AE1A4A95FF82413), // 4.560662751518062e-25
        ],
    },
];

impl crate::generic::Erf for f64 {
    #[inline]
    fn erf_finite(x: Self) -> Self {
        let ax = x.abs();
        if ax < POLY_MAX {
            if x.exponent() < -960 {
                // erf(x) = x * 2 / sqrt(π) - ..., where the next terms are
                // negligible, and the result can be subnormal
                mul_tiny(x, FRAC_2_SQRT_PI)
            } else {
                let (hi, lo) = erf_small(x);
                hi + lo
            }
        } else if ax < PLAIN_MIN {
            // erf(x) = (1 - erfc(|x|)) * sgn(x), where erfc(|x|) < 0.48
            let (hi, lo, _) = erfc_large(ax);
            let s = F64x2::fast_add11(1.0, -hi);
            (s.hi() + (s.lo() - lo)).copysign(x)
        } else if ax < SATURATION_MIN {
            // erfc(|x|) < 2^-20, see `PLAIN_MIN`
            (1.0 - erfc_plain(ax)).copysign(x)
        } else {
            1.0f64.copysign(x)
        }
    }

    #[inline]
    fn erfc_finite(x: Self) -> Self {
        let ax = x.abs();
        if ax < POLY_MAX {
            if x.exponent() < -55 {
                // |erf(x)| < 2^-54, so erfc(x) = 1 - erf(x) rounds to 1
                return 1.0;
            }
            // erfc(x) = 1 - erf(x), where |erf(x)| < 0.53
            let (hi, lo) = erf_small(x);
            let s = F64x2::fast_add11(1.0, -hi);
            s.hi() + (s.lo() - lo)
        } else if x.is_sign_negative() {
            if ax < PLAIN_MIN {
                // erfc(x) = 2 - erfc(|x|), where erfc(|x|) < 0.48
                let (hi, lo, _) = erfc_large(ax);
                let s = F64x2::fast_add11(2.0, -hi);
                s.hi() + (s.lo() - lo)
            } else if ax < SATURATION_MIN {
                // erfc(|x|) < 2^-20, see `PLAIN_MIN`
                2.0 - erfc_plain(ax)
            } else {
                2.0
            }
        } else if x < ERFC_ZERO_MIN {
            let (hi, lo, edelta) = erfc_large(x);
            if edelta == 0 {
                hi + lo
            } else {
                // The result can be subnormal, scale it rounding only once
                F64x2::fast_add11(hi, lo).scalbn_to_f64(edelta)
            }
        } else {
            0.0
        }
    }

    #[inline]
    fn erfcx_finite(x: Self) -> Self {
        let ax = x.abs();
        if ax < POLY_MAX {
            if x.exponent() < -55 {
                // |erf(x)| < 2^-54 and x^2 < 2^-110, so
                // erfcx(x) = exp(x^2) * (1 - erf(x)) rounds to 1
                1.0
            } else {
                erfcx_small(x)
            }
        } else if x.is_sign_negative() {
            if ax < ERFCX_INF_MIN {
                erfcx_neg(ax)
            } else {
                f64::INFINITY
            }
        } else if x < ASYMPT_MIN {
            let (hi, lo) = erfcx_poly(x);
            hi + lo
        } else {
            erfcx_asympt(x)
        }
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= erf(x)`, for
/// `2^-960 <= |x| < POLY_MAX`, with the relative error of `erf_poly`.
///
/// `|lo| <= 2^-22 * |hi|`.
#[inline]
fn erf_small(x: f64) -> (f64, f64) {
    // erf(x) = x * P(x^2)
    let (p_hi, p_lo) = if x.exponent() < -35 {
        // P(x^2) = 2 / sqrt(π) * (1 - x^2 / 3 + ...), where the second term
        // is less than 2^-71.5
        (FRAC_2_SQRT_PI.hi(), FRAC_2_SQRT_PI.lo())
    } else {
        // u = x^2 as a normalized double-word
        let (u_hi, u_lo) = square_parts(x);
        erf_poly(u_hi, u_lo)
    };
    Split::new(x).mul(p_hi, p_lo)
}

/// Returns `(hi, lo)` such that `hi + lo ~= P(u) = erf(x) / x`, where
/// `u = x^2 = u_hi + u_lo` (with `|u_lo| <= 2^-53 * u_hi`), for
/// `0 <= u < 1/4`, with a relative error of about 2^-68.5 (mostly from the
/// evaluation of the terms of degree 4 and higher).
///
/// `|lo| <= 2^-23 * |hi|`.
#[inline]
fn erf_poly(u_hi: f64, u_lo: f64) -> (f64, f64) {
    // P(u) ~= K0 + K1 * u + ... + K10 * u^10, with a relative error of
    // 2^-72.8
    // GENERATE: erf_poly F64x2:4,f64 11 0.5
    const K0: F64x2 = F64x2::from_bits(0x3FF20DD750429B6D, 0x3C71AE30C3E6BAB7); // 1.128379167095512573896029169617e0
    const K1: F64x2 = F64x2::from_bits(0xBFD812746B0379E7, 0x3C6F2AEA2CA8AC65); // -3.761263890318375245071428308365e-1
    const K2: F64x2 = F64x2::from_bits(0x3FBCE2F21A042BE0, 0x3C51AEDF43F3BD05); // 1.128379167095512374860192791321e-1
    const K3: F64x2 = F64x2::from_bits(0xBF9B82CE312889EC, 0x3BF44B4B6DD7D974); // -2.686617064513001355904016258331e-2
    const K4: f64 = f64::from_bits(0x3F7565BCD0E5F32C); // 5.2239776254026475e-3
    const K5: f64 = f64::from_bits(0xBF4C02DB3F9C4818); // -8.54832701607802e-4
    const K6: f64 = f64::from_bits(0x3F1F9A3249DFA211); // 1.2055332123817759e-4
    const K7: f64 = f64::from_bits(0xBEEF4D1CF0A3CC89); // -1.4925586146996457e-5
    const K8: f64 = f64::from_bits(0x3EBB9D186BCF2B75); // 1.6459020308774685e-6
    const K9: f64 = f64::from_bits(0xBE85D75140791251); // -1.6272875334918868e-7
    const K10: f64 = f64::from_bits(0x3E4C5F86436CD7A2); // 1.3212274036643748e-8

    // The terms of degree 4 and higher (less than 2^-15.5 relative to
    // `P(u)`) are evaluated with plain `f64` arithmetic, where the last step
    // is in Horner form, so the relative error of `v` is about 2^-52.8
    let u = u_hi;
    let u2 = u * u;
    let u4 = u2 * u2;
    let t = ((K5 + u * K6) + u2 * (K7 + u * K8)) + u4 * (K9 + u * K10);
    let v = K4 + u * t;

    // The next steps use exact products
    let us = Split::new(u);
    let (h, l) = us.mul_add(v, 0.0, K3);
    let (h, l) = us.mul_add(h, l, K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (p_hi, p_lo) = us.mul_add(h, l, K0);

    // P(u_hi + u_lo) ~= P(u_hi) + u_lo * P'(u_hi), where
    // P'(u) ~= K1 + 2 * K2 * u + 3 * K3 * u^2 + 4 * K4 * u^3 (the next
    // terms are less than 2^-14 relative to `P'(u)`, and
    // |u_lo * P'(u)| < 2^-56.5 * P(u))
    let dp = K1.hi() + u * (2.0 * K2.hi() + u * (3.0 * K3.hi() + u * (4.0 * K4)));
    (p_hi, p_lo + u_lo * dp)
}

/// Returns `(hi, lo, edelta)` such that `(hi + lo) * 2^edelta ~= erfc(x)`,
/// for `POLY_MAX <= x < 28`, with a relative error of about 2^-67:
/// * `exp(-x^2)`: 2^-68 (see `Reduced::eval`).
/// * `erfcx(x)`: 2^-67.5 (see `erfcx_poly`).
///
/// `edelta` is the one of `Reduced::eval_scaled` for `exp(-x^2)`, which is
/// zero when `x < 25.79`:
/// * When it is zero, the result is greater than 2^-966, so the rounding of
///   the terms of its products that are subnormal (at most 2^-1075) is
///   negligible.
/// * When it is not zero, `2^-6 < hi < 2^-4` and the result can be subnormal.
///
/// `|lo| <= 2^-16 * |hi|`.
#[inline]
fn erfc_large(x: f64) -> (f64, f64, i32) {
    // x^2 = p + pe, where the error of `square_parts` is negligible
    let (p, pe) = square_parts(x);

    // exp(-x^2) = (e_hi + e_lo) * 2^edelta
    let red = Reduced::exp_sum(-p, -pe);
    let (e_hi, e_lo, edelta) = red.eval_scaled(exp_poly(red.r));

    // erfcx(x) = c_hi + c_lo
    let (c_hi, c_lo) = erfcx_poly(x);

    // erfc(x) = e_hi * (c_hi + c_lo) + e_lo * (c_hi + c_lo), where
    // `|e_lo| < 2^-17 * e_hi`, so the rounding errors of the second product
    // are less than 2^-69 relative to the result
    let (hi, lo) = Split::new(e_hi).mul(c_hi, c_lo);
    (hi, lo + e_lo * (c_hi + c_lo), edelta)
}

/// Returns `erfc(x)`, for `POLY_MAX <= x < 25` (where `exp(-x^2)` does not
/// need to be scaled, see `erfc_large`), with a relative error of about
/// 2^-50.5 (mostly from `erfcx_plain`).
///
/// It is much faster than `erfc_large`, for the results that are subtracted
/// from much larger values.
#[inline]
fn erfc_plain(x: f64) -> f64 {
    // exp(-x^2) = e_hi + e_lo, like in `erfc_large`
    let (p, pe) = square_parts(x);
    let red = Reduced::exp_sum(-p, -pe);
    let (e_hi, e_lo) = red.eval(exp_poly(red.r));

    (e_hi + e_lo) * erfcx_plain(x)
}

/// Returns `erfcx(x) = exp(x^2) * (1 - erf(x))` rounded to `f64`, for
/// `2^-960 <= |x| < POLY_MAX`.
///
/// The relative error before the final rounding is about 2^-67.5 (see
/// `erf_small` and `Reduced::eval`).
#[inline]
fn erfcx_small(x: f64) -> f64 {
    // c_hi + c_lo = 1 - erf(x), where |erf(x)| < 0.53
    let (hi, lo) = erf_small(x);
    let s = F64x2::fast_add11(1.0, -hi);
    let (c_hi, c_lo) = (s.hi(), s.lo() - lo);

    // e_hi + e_lo = exp(x^2)
    let (p, pe) = square_parts(x);
    let red = Reduced::exp_sum(p, pe);
    let (e_hi, e_lo) = red.eval(exp_poly(red.r));

    // The product is calculated like in `erfc_large`
    let (hi, lo) = Split::new(e_hi).mul(c_hi, c_lo);
    hi + (lo + e_lo * (c_hi + c_lo))
}

/// Returns `erfcx(-x) = 2 * exp(x^2) - erfcx(x)` rounded to `f64`, for
/// `POLY_MAX <= x < 32`, where `erfcx(x) < 0.62` and `2 * exp(x^2) > 2.56`.
///
/// The relative error before the final rounding is about 2^-68 (see
/// `Reduced::eval` and `erfcx_poly`). The result overflows when `x > ~26.63`.
#[inline]
fn erfcx_neg(x: f64) -> f64 {
    // x^2 = p + pe, where the error of `square_parts` is negligible
    let (p, pe) = square_parts(x);
    let red = Reduced::exp_sum(p, pe).mul_exp2(1);
    if x < ERFCX_EXP_MIN {
        // 2 * exp(x^2) = e_hi + e_lo, where `e_hi` is exact
        let (e_hi, e_lo) = red.eval(exp_poly(red.r));
        if x < PLAIN_MIN {
            let (c_hi, c_lo) = erfcx_poly(x);
            let s = F64x2::fast_add11(e_hi, -c_hi);
            s.hi() + (s.lo() + (e_lo - c_lo))
        } else {
            // erfcx(x) < 2^-20 * exp(x^2), see `PLAIN_MIN`
            e_hi + (e_lo - erfcx_plain(x))
        }
    } else {
        // erfcx(x) is negligible, see `ERFCX_EXP_MIN`
        red.eval_to_f64(exp_poly(red.r))
    }
}

/// Returns `erfcx(x)` rounded to `f64`, for `x >= ASYMPT_MIN`.
///
/// The relative error before the final rounding is about 2^-70 (mostly from
/// the evaluation of the terms of degree 2 and higher of the polynomial).
#[inline]
fn erfcx_asympt(x: f64) -> f64 {
    // erfcx(x) ~= Q(s) / x, with s = 1 / x^2 and
    // Q(s) = K0 + K1 * s + ... + K7 * s^7, with a relative error of 2^-78.9
    // GENERATE: erfcx_asympt_poly F64x2:2,f64 8 28
    const K0: F64x2 = F64x2::from_bits(0x3FE20DD750429B6D, 0x3C61AE3A914FED80); // 5.641895835477562869480794515608e-1
    const K1: F64x2 = F64x2::from_bits(0xBFD20DD750429B6D, 0xBC51642D638B269B); // -2.820947917738781434113156827457e-1
    const K2: f64 = f64::from_bits(0x3FDB14C2F863E90D); // 4.2314218766081596e-1
    const K3: f64 = f64::from_bits(0xBFF0ECF9DB3DDBDC); // -1.057855469143525e0
    const K4: f64 = f64::from_bits(0x400D9EB53C005CA9); // 3.702494114647909e0
    const K5: f64 = f64::from_bits(0xC030A942DD4A3130); // -1.6661176519978824e1
    const K6: f64 = f64::from_bits(0x4056E5E4F81D2600); // 9.159210017057922e1
    const K7: f64 = f64::from_bits(0xC081ED1B9946B3CF); // -5.736384759448646e2

    if x.exponent() >= 35 {
        // Q(s) = K0 * (1 - s / 2 + ...), where the terms after the first are
        // less than 2^-70. `x` is scaled so its reciprocal is not subnormal,
        // and the result, which is subnormal when `x > ~2^1021.2`, is rounded
        // only once.
        let (t_hi, t_lo) = recip_parts(x * f64::exp2i_fast(-512));
        return (F64x2::fast_add11(t_hi, t_lo) * K0).scalbn_to_f64(-512);
    }

    // 1 / x = t_hi + t_lo
    let (t_hi, t_lo) = recip_parts(x);
    // s = s_hi + s_lo, where `s_hi` is `t_hi^2` rounded
    let (s, se) = square_parts(t_hi);
    let s_lo = se + 2.0 * t_hi * t_lo;

    // The terms of degree 2 and higher (less than 2^-19.5 relative to
    // `Q(s)`) are evaluated with plain `f64` arithmetic
    let s2 = s * s;
    let s4 = s2 * s2;
    let w = ((K2 + s * K3) + s2 * (K4 + s * K5)) + s4 * (K6 + s * K7);

    // The next steps use exact products
    let ss = Split::new(s);
    let (h, l) = ss.mul_add(w, 0.0, K1);
    let (q_hi, q_lo) = ss.mul_add(h, l, K0);
    // Q(s_hi + s_lo) ~= Q(s_hi) + s_lo * K1 (the next terms are negligible)
    let q_lo = q_lo + s_lo * K1.hi();

    // The product is calculated like in `erfc_large`
    let (hi, lo) = Split::new(t_hi).mul(q_hi, q_lo);
    hi + (lo + t_lo * (q_hi + q_lo))
}

/// Returns `(d, poly)`, where `poly` is the approximation of `erfcx` in the
/// interval of `x` and `d = x - c` (see `ErfcxPoly`), for
/// `POLY_MAX <= x < 28`.
///
/// `|d| <= c / 9`.
#[inline]
fn erfcx_interval(x: f64) -> (f64, &'static ErfcxPoly) {
    // The interval of `x` is `[a, a + w)`, where `a` is `x` with all but the
    // two leading bits of the mantissa cleared and `w` is a quarter of the
    // binade, so its center is `a` with the next bit set
    let bits = x.to_bits();
    let i = (bits >> ERFCX_TBL_SHIFT) - ERFCX_TBL_FIRST;
    let c = f64::from_bits((bits & (u64::MAX << ERFCX_TBL_SHIFT)) | (1 << (ERFCX_TBL_SHIFT - 1)));
    // exact
    (x - c, &ERFCX_TBL[i as usize])
}

/// Returns `(hi, lo)` such that `hi + lo ~= erfcx(x)`, for
/// `POLY_MAX <= x < 28`, with a relative error of about 2^-67.5:
/// * Polynomial approximation: less than 2^-68.5 (2^-70 without the rounding
///   of the coefficients).
/// * Evaluation of the terms of degree 5 and higher (less than 2^-15.8
///   relative to the result): ~2^-68.5.
///
/// `|lo| <= 2^-23 * |hi|`.
#[inline]
fn erfcx_poly(x: f64) -> (f64, f64) {
    let (d, poly) = erfcx_interval(x);
    let [k0, k1, k2, k3, k4] = poly.k_hi;
    let [k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16] = poly.k_lo;

    // The terms of degree 5 and higher are evaluated with plain `f64`
    // arithmetic, where the last step is in Horner form, so the relative
    // error of `v` is about 2^-52.8
    let d2 = d * d;
    let d4 = d2 * d2;
    let d8 = d4 * d4;
    let t = (((k6 + d * k7) + d2 * (k8 + d * k9)) + d4 * ((k10 + d * k11) + d2 * (k12 + d * k13)))
        + d8 * ((k14 + d * k15) + d2 * k16);
    let v = k5 + d * t;

    // The next steps use exact products
    let ds = Split::new(d);
    let (h, l) = ds.mul_add(v, 0.0, k4);
    let (h, l) = ds.mul_add(h, l, k3);
    let (h, l) = ds.mul_add(h, l, k2);
    let (h, l) = ds.mul_add(h, l, k1);
    ds.mul_add(h, l, k0)
}

/// Returns `erfcx(x)`, for `POLY_MAX <= x < 28`, with a relative error of
/// about 2^-51.
///
/// It evaluates the whole polynomial of `ERFCX_TBL` with plain `f64`
/// arithmetic.
#[inline]
fn erfcx_plain(x: f64) -> f64 {
    let (d, poly) = erfcx_interval(x);
    let [k0, k1, k2, k3, k4] = poly.k_hi.map(F64x2::hi);
    let [k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16] = poly.k_lo;

    let d2 = d * d;
    let d4 = d2 * d2;
    let d8 = d4 * d4;
    let t = (((k8 + d * k9) + d2 * (k10 + d * k11))
        + d4 * ((k12 + d * k13) + d2 * (k14 + d * k15)))
        + d8 * k16;
    (((k0 + d * k1) + d2 * (k2 + d * k3)) + d4 * ((k4 + d * k5) + d2 * (k6 + d * k7))) + d8 * t
}

/// Returns `erfcx(x)`, for `POLY_MAX <= x < 28`, with a relative error of
/// about 2^-40.5.
///
/// It is meant for the `f32` functions: it uses plain `f64` arithmetic, and
/// only the terms up to degree 12 of the polynomials of `ERFCX_TBL` (the next
/// ones are less than 2^-40.5 relative to the result).
#[inline]
pub(crate) fn erfcx_short(x: f64) -> f64 {
    let (d, poly) = erfcx_interval(x);
    let [k0, k1, k2, k3, k4] = poly.k_hi.map(F64x2::hi);
    let [k5, k6, k7, k8, k9, k10, k11, k12, ..] = poly.k_lo;

    let d2 = d * d;
    let d4 = d2 * d2;
    let d8 = d4 * d4;
    (((k0 + d * k1) + d2 * (k2 + d * k3)) + d4 * ((k4 + d * k5) + d2 * (k6 + d * k7)))
        + d8 * (((k8 + d * k9) + d2 * (k10 + d * k11)) + d4 * k12)
}

#[cfg(test)]
mod tests {
    use super::{ASYMPT_MIN, ERFCX_TBL, ERFCX_TBL_FIRST, ERFCX_TBL_SHIFT, erfcx_interval};

    const PREC: u32 = 256;

    /// Returns `log2` of the relative error of `approx` with respect to
    /// `exact`.
    fn rel_error_log2(approx: &rug::Float, exact: &rug::Float) -> f64 {
        let err = rug::Float::with_val(PREC, approx - exact) / exact;
        err.abs().log2().to_f64()
    }

    #[test]
    fn test_erfcx_tbl() {
        // `erfcx_poly` requires a relative error less than 2^-68.5 from the
        // polynomials. With only the high part of the coefficients,
        // `erfcx_plain` requires a relative error less than 2^-52, and
        // `erfcx_short` requires a relative error less than 2^-40.5 from
        // their terms up to degree 12.
        let limit = |i: usize| f64::from_bits((i as u64 + ERFCX_TBL_FIRST) << ERFCX_TBL_SHIFT);
        assert_eq!(limit(0), 0.5);
        assert_eq!(limit(ERFCX_TBL.len()), ASYMPT_MIN);

        for (i, poly) in ERFCX_TBL.iter().enumerate() {
            // The interval is `[a, b)`
            let (a, b) = (limit(i), limit(i + 1));
            let below_b = f64::from_bits(b.to_bits() - 1);
            let xs = (0..64)
                .map(|j| a + (b - a) * (f64::from(j) / 64.0))
                .chain([below_b]);
            for x in xs {
                let (d, x_poly) = erfcx_interval(x);
                assert!(core::ptr::eq(x_poly, poly));
                assert_eq!(d, x - (a + b) / 2.0);

                // erfcx(x) = exp(x^2) * erfc(x)
                let xr = rug::Float::with_val(PREC, x);
                let exact = rug::Float::with_val(PREC, xr.square_ref()).exp() * xr.erfc();

                // Evaluates the first `num` terms of the polynomial, with or
                // without the low parts of the double-word coefficients
                let eval = |num: usize, with_lo: bool| {
                    let k_hi = poly.k_hi.iter().map(|k| {
                        let k_lo = if with_lo { k.lo() } else { 0.0 };
                        rug::Float::with_val(PREC, k.hi()) + k_lo
                    });
                    let k_lo = poly.k_lo.iter().map(|&k| rug::Float::with_val(PREC, k));
                    let skip = poly.k_hi.len() + poly.k_lo.len() - num;
                    (k_hi.chain(k_lo).rev().skip(skip))
                        .fold(rug::Float::new(PREC), |acc, k| acc * d + k)
                };
                assert!(rel_error_log2(&eval(17, true), &exact) < -68.5, "x = {x:e}");
                assert!(
                    rel_error_log2(&eval(17, false), &exact) < -52.0,
                    "x = {x:e}"
                );
                assert!(
                    rel_error_log2(&eval(13, false), &exact) < -40.5,
                    "x = {x:e}"
                );
            }
        }
    }
}
