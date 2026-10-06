//! Exponential functions for `f64`, and the parts that `f32` shares.
//!
//! The argument is reduced to `x = m * ln(2) / N + r` (or the equivalent for
//! the other bases), where `N = 2^EXP2_TBL_BITS`, `m = k * N + j` is an
//! integer with `0 <= j < N`, and `|r| <= ~ln(2) / (2 * N)`, so
//!
//! `exp(x) = 2^k * 2^(j / N) * exp(r)`
//!
//! `T`, `2^(j / N)` rounded to `f64`, is taken from `EXP2_TBL`, and
//! `S = 2^k * T` is obtained by adding `k` to the exponent of `T` (see
//! `exp2_tbl`). `exp(r) - 1 - r` is approximated with a short polynomial
//! `r^2 * P(r)`.
//!
//! To keep the error close to 0.5 ULP, the leading part of the result is
//! calculated exactly: `2^(m / N) = Sh + U`, where `Sh` is `S` truncated to 25
//! significant bits and `U = S * υ`, with `υ` also taken from `EXP2_TBL`, and
//! `r = a + b`, where `a` is a multiple of 2^-27. Then `1 + a` has at most 28
//! bits, `Sh * (1 + a)` is exact, and
//!
//! `2^(m / N) * exp(r) ~= Sh * (1 + a) + (Sh * b + U * (1 + r) + S * r^2 * P(r))`
//!
//! The second part is small, so its rounding errors are small relative to the
//! result, and the result is rounded only once, when both parts are added.
//!
//! Reductions:
//! * `exp`: Cody-Waite with a split of `ln(2) / N`.
//! * `exp2`: `x * N = m + f`, with exact `f`, and `r = f * ln(2) / N`.
//! * `exp10`: Cody-Waite with a split of `log10(2) / N`, giving `t`, and
//!   `r = t * ln(10)`.
//!
//! In the latter two, `r` is calculated as a sum of two `f64`, where the
//! product with the high part of `ln(2) / N` or `ln(10)` is exact.

use super::f64x2::F64x2;
use super::round_i32;
use crate::generic::scalbn_medium;
use crate::traits::Float as _;

/// Number of bits of the index of `EXP2_TBL`
pub(crate) const EXP2_TBL_BITS: u32 = {
    assert!(EXP2_TBL.len().is_power_of_two());
    EXP2_TBL.len().ilog2()
};

/// `N`, the number of entries of `EXP2_TBL`
const TBL_N: f64 = (1 << EXP2_TBL_BITS) as f64;

const TBL_MASK: i32 = (1 << EXP2_TBL_BITS) - 1;

// GENERATE: consts f64 LOG2_E LOG2_10 LN_2 LN_10
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LOG2_10: f64 = f64::from_bits(0x400A934F0979A371); // 3.321928094887362e0
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1
const LN_10: f64 = f64::from_bits(0x40026BB1BBB55516); // 2.302585092994046e0

// `ln(2) / N` and `log10(2) / N`, the first parts have 35 bits, so their
// products with integers of up to 18 bits are exact
// GENERATE: split_const LN_2_N LN_2 -7 35 53
const LN_2_N0: f64 = f64::from_bits(0x3F762E42FEF80000); // 5.415212347998022e-3
const LN_2_N1: f64 = f64::from_bits(0x3D41CF79ABC9E3B4); // 1.2655086083325438e-13

// GENERATE: split_const LOG10_2_N LOG10_2 -7 35 53
const LOG10_2_N0: f64 = f64::from_bits(0x3F634413509C0000); // 2.3517968410260437e-3
const LOG10_2_N1: f64 = f64::from_bits(0x3D3BCFF7988F895A); // 9.880939360664394e-14

// `ln(2)` and `ln(10)`, the first parts have 27 bits, so their products with
// numbers of up to 26 bits are exact
// GENERATE: split_const LN_2_X LN_2 0 27 53
const LN_2_X0: f64 = f64::from_bits(0x3FE62E42FC000000); // 6.93147175014019e-1
const LN_2_X1: f64 = f64::from_bits(0x3E37D1CF79ABC9E4); // 5.5459262969660605e-9

// GENERATE: split_const LN_10_X LN_10 0 27 53
const LN_10_X0: f64 = f64::from_bits(0x40026BB1B8000000); // 2.3025850653648376e0
const LN_10_X1: f64 = f64::from_bits(0x3E5DAAA8AC16EA57); // 2.7629208037533617e-8

// GENERATE: exp2_table 7 25
// EXP2_TBL[i] = (bits(T), bits(υ)), where T is 2^(i / 128) rounded,
// Th is T truncated to 25 bits and 2^(i / 128) = Th + T * υ
static EXP2_TBL: [(u64, u64); 128] = [
    (0x3FF0000000000000, 0x0000000000000000), // 1e0
    (0x3FF0163DA9FB3335, 0x3E63DACD0BCE42E6), // 1.0054299011128027e0
    (0x3FF02C9A3E778061, 0x3E6C9F36FE711C01), // 1.0108892860517005e0
    (0x3FF04315E86E7F85, 0x3E60976E34393F93), // 1.016378314910953e0
    (0x3FF059B0D3158574, 0x3E4824D3F33CA46E), // 1.0218971486541166e0
    (0x3FF0706B29DDF6DE, 0x3E6334FA3F470AAF), // 1.0274459491187637e0
    (0x3FF0874518759BC8, 0x3E6060C07DA5BEB8), // 1.0330248790212284e0
    (0x3FF09E3ECAC6F383, 0x3E64C0A7228B3538), // 1.0386341019613787e0
    (0x3FF0B5586CF9890F, 0x3E68D96D26E5A82E), // 1.0442737824274138e0
    (0x3FF0CC922B7247F7, 0x3E65CDC8FB01E3E9), // 1.0499440858006872e0
    (0x3FF0E3EC32D3D1A2, 0x3E456D51C74D1E55), // 1.0556451783605572e0
    (0x3FF0FB66AFFED31B, 0x3E6E240F17B5C9EA), // 1.061377227289262e0
    (0x3FF11301D0125B51, 0x3DF133A6ECF5608F), // 1.0671404006768237e0
    (0x3FF12ABDC06C31CC, 0x3E1935BEFF7E832A), // 1.0729348675259756e0
    (0x3FF1429AAEA92DE0, 0x3E6B2E50DD7BBD6B), // 1.0787607977571199e0
    (0x3FF15A98C8A58E51, 0x3E5FE3723D7DEEB1), // 1.0846183622133092e0
    (0x3FF172B83C7D517B, 0x3E66E7E8574DDAE0), // 1.0905077326652577e0
    (0x3FF18AF9388C8DEA, 0x3E5F304C7D7B92F7), // 1.0964290818163769e0
    (0x3FF1A35BEB6FCB75, 0x3E64BFC20CD7FF6D), // 1.102382583307841e0
    (0x3FF1BBE084045CD4, 0x3E4CFE886209F247), // 1.1083684117236787e0
    (0x3FF1D4873168B9AA, 0x3E343B2CD6416D45), // 1.1143867425958924e0
    (0x3FF1ED5022FCD91D, 0x3E455510F37B279B), // 1.1204377524096067e0
    (0x3FF2063B88628CD6, 0x3E5DC5DE8C7C1DB2), // 1.1265216186082418e0
    (0x3FF21F49917DDC96, 0x3E35124B99949265), // 1.1326385195987192e0
    (0x3FF2387A6E756238, 0x3E69649025814D7C), // 1.1387886347566916e0
    (0x3FF251CE4FB2A63F, 0x3E6B6BA54D6133DC), // 1.1449721444318042e0
    (0x3FF26B4565E27CDD, 0x3E54728B55DE778C), // 1.1511892299529827e0
    (0x3FF284DFE1F56381, 0x3E3B130062B7BD3A), // 1.1574400736337511e0
    (0x3FF29E9DF51FDEE1, 0x3E519D365BE83037), // 1.1637248587775775e0
    (0x3FF2B87FD0DAD990, 0x3E276169C105C1FE), // 1.1700437696832502e0
    (0x3FF2D285A6E4030B, 0x3E576E04015DDF5E), // 1.1763969916502812e0
    (0x3FF2ECAFA93E2F56, 0x3E5F4214DBBE81BC), // 1.182784710984341e0
    (0x3FF306FE0A31B715, 0x3E61250010DE5BC2), // 1.189207115002721e0
    (0x3FF32170FC4CD831, 0x3E64931AB02C4510), // 1.1956643920398273e0
    (0x3FF33C08B26416FF, 0x3E3FD28E60AA99CB), // 1.202156731452703e0
    (0x3FF356C55F929FF1, 0x3E69C4A29B6CF23B), // 1.2086843236265816e0
    (0x3FF371A7373AA9CB, 0x3E57CB7A75B9478E), // 1.215247359980469e0
    (0x3FF38CAE6D05D866, 0x3E655110171CB412), // 1.2218460329727576e0
    (0x3FF3A7DB34E59FF7, 0x3E4FE3BE1B4D75C6), // 1.22848053610687e0
    (0x3FF3C32DC313A8E5, 0x3E43EDA0127C96D8), // 1.2351510639369334e0
    (0x3FF3DEA64C123422, 0x3E6370BE3B9F07C8), // 1.241857812073484e0
    (0x3FF3FA4504AC801C, 0x3E4DF22ECDEE5E4D), // 1.2486009771892048e0
    (0x3FF4160A21F72E2A, 0x3E390D1A31B45BC8), // 1.255380757024691e0
    (0x3FF431F5D950A897, 0x3E5D852A5AF3070D), // 1.2621973503942507e0
    (0x3FF44E086061892D, 0x3E1336DE2B9BDFC9), // 1.2690509571917332e0
    (0x3FF46A41ED1D0057, 0x3E648DFF677D77C9), // 1.275941778396392e0
    (0x3FF486A2B5C13CD0, 0x3E51F191514161CF), // 1.2828700160787783e0
    (0x3FF4A32AF0D7D3DE, 0x3E24EA8B65082ACA), // 1.2898358734066657e0
    (0x3FF4BFDAD5362A27, 0x3E50132059320084), // 1.2968395546510096e0
    (0x3FF4DCB299FDDD0D, 0x3E5EA6EBE7DFA689), // 1.3038812651919358e0
    (0x3FF4F9B2769D2CA7, 0x3E542E35203782D7), // 1.3109612115247644e0
    (0x3FF516DAA2CF6642, 0x3E410E5956F4E0BB), // 1.318079601266064e0
    (0x3FF5342B569D4F82, 0x3E53F6F7DA076891), // 1.3252366431597413e0
    (0x3FF551A4CA5D920F, 0x3E5F1E17EAA5EB79), // 1.3324325470831615e0
    (0x3FF56F4736B527DA, 0x3E54071C2650EC31), // 1.339667524053303e0
    (0x3FF58D12D497C7FD, 0x3E4B47698373412E), // 1.3469417862329458e0
    (0x3FF5AB07DD485429, 0x3E639DB1C4F3D6CF), // 1.3542555469368927e0
    (0x3FF5C9268A5946B7, 0x3E5E66C7F4A6941B), // 1.3616090206382248e0
    (0x3FF5E76F15AD2148, 0x3E5095CE9660151F), // 1.3690024229745905e0
    (0x3FF605E1B976DC09, 0x3E5B80F723BCF3DD), // 1.3764359707545302e0
    (0x3FF6247EB03A5585, 0x3E05136894DCD33F), // 1.383909881963832e0
    (0x3FF6434634CCC320, 0x3E4B98C5EFF0C7A5), // 1.3914243757719262e0
    (0x3FF6623882552225, 0x3E3AAD5BD0788FBD), // 1.3989796725383112e0
    (0x3FF68155D44CA973, 0x3E4874179507C12A), // 1.4065759938190154e0
    (0x3FF6A09E667F3BCD, 0x3E526055C2A3609F), // 1.4142135623730951e0
    (0x3FF6C012750BDABF, 0x3E4C645DB064D097), // 1.4218926021691656e0
    (0x3FF6DFB23C651A2F, 0x3E61571B07BF84DC), // 1.42961333839197e0
    (0x3FF6FF7DF9519484, 0x3E59EEB56D9EEA04), // 1.4373759974489824e0
    (0x3FF71F75E8EC5F74, 0x3E58B2BB7B9A2A7D), // 1.4451808069770467e0
    (0x3FF73F9A48A58174, 0x3E57CD7F667DB4CA), // 1.4530279958490526e0
    (0x3FF75FEB564267C9, 0x3E512364CE10E1CD), // 1.460917794180647e0
    (0x3FF780694FDE5D3F, 0x3E659B5A0FF03B44), // 1.4688504333369818e0
    (0x3FF7A11473EB0187, 0x3E45394E39334248), // 1.4768261459394993e0
    (0x3FF7C1ED0130C132, 0x3E29A7C9F7AAB01B), // 1.4848451658727524e0
    (0x3FF7E2F336CF4E62, 0x3E523EE6BF64C2D3), // 1.4929077282912648e0
    (0x3FF80427543E1A12, 0x3E469CA135B2A902), // 1.5010140696264256e0
    (0x3FF82589994CCE13, 0x3E58A6428DA911B6), // 1.5091644275934228e0
    (0x3FF8471A4623C7AD, 0x3E502F7656ADADD7), // 1.5173590411982147e0
    (0x3FF868D99B4492ED, 0x3E5D8B2109194C01), // 1.5255981507445384e0
    (0x3FF88AC7D98A6699, 0x3E58E133840ACD70), // 1.533881997840956e0
    (0x3FF8ACE5422AA0DB, 0x3E367A1CA0DD3CB5), // 1.5422108254079407e0
    (0x3FF8CF3216B5448C, 0x3E514DFCB413F3E0), // 1.550584877685e0
    (0x3FF8F1AE99157736, 0x3E574E8A5547739D), // 1.559004400237837e0
    (0x3FF9145B0B91FFC6, 0x3E5D86AD43FC41CC), // 1.567469639965553e0
    (0x3FF93737B0CDC5E5, 0x3E20522EEA80765A), // 1.5759808451078865e0
    (0x3FF95A44CBC8520F, 0x3E5DBE6232C520C5), // 1.5845382652524937e0
    (0x3FF97D829FDE4E50, 0x3E63EBBD9BBC7F0D), // 1.593142151342267e0
    (0x3FF9A0F170CA07BA, 0x3E1F882997ED8FF9), // 1.6017927556826934e0
    (0x3FF9C49182A3F090, 0x3E3A3B5E32F76D6B), // 1.6104903319492543e0
    (0x3FF9E86319E32323, 0x3E586CAC7A8F2DF1), // 1.6192351351948637e0
    (0x3FFA0C667B5DE565, 0x3E5BED80594E7ADF), // 1.6280274218573478e0
    (0x3FFA309BEC4A2D33, 0x3E5E084965C81C19), // 1.6368674497669644e0
    (0x3FFA5503B23E255D, 0x3E35CDD5B1362D2A), // 1.645755478153965e0
    (0x3FFA799E1330B358, 0x3E3ED910BD6F9180), // 1.6546917676561943e0
    (0x3FFA9E6B5579FDBF, 0x3E4A55A712F0BB3F), // 1.6636765803267364e0
    (0x3FFAC36BBFD3F37A, 0x3E62ECC60BFEB62D), // 1.6727101796415966e0
    (0x3FFAE89F995AD3AD, 0x3E563FE3EA2EF37C), // 1.681792830507429e0
    (0x3FFB0E07298DB666, 0x3E569980AF4B2E43), // 1.6909247992693053e0
    (0x3FFB33A2B84F15FB, 0x3E538C9845CA423B), // 1.7001063537185235e0
    (0x3FFB59728DE5593A, 0x3E604240E4ADB5B3), // 1.709337763100463e0
    (0x3FFB7F76F2FB5E47, 0x3E3BC2CE18032F5B), // 1.718619298122478e0
    (0x3FFBA5B030A1064A, 0x3E174C0A7797CADF), // 1.7279512309618377e0
    (0x3FFBCC1E904BC1D2, 0x3E05CD7D4B284B28), // 1.7373338352737062e0
    (0x3FFBF2C25BD71E09, 0x3E5B1D16FC795902), // 1.746767386199169e0
    (0x3FFC199BDD85529C, 0x3E5ECB6E0D6D1973), // 1.7562521603732995e0
    (0x3FFC40AB5FFFD07A, 0x3E621F13BE4B45DC), // 1.7657884359332727e0
    (0x3FFC67F12E57D14B, 0x3E602861C14631EF), // 1.7753764925265212e0
    (0x3FFC8F6D9406E7B5, 0x3E420C428F86D9D8), // 1.785016611318935e0
    (0x3FFCB720DCEF9069, 0x3E5CD4B9ECB6ED49), // 1.7947090750031072e0
    (0x3FFCDF0B555DC3FA, 0x3E47CA8DA359ABFE), // 1.804454167806624e0
    (0x3FFD072D4A07897C, 0x3E561CD0F4A685D6), // 1.8142521755003989e0
    (0x3FFD2F87080D89F2, 0x3E51A8A98B1F72C2), // 1.8241033854070534e0
    (0x3FFD5818DCFBA487, 0x3E5C50EA85A92131), // 1.8340080864093424e0
    (0x3FFD80E316C98398, 0x3E4D722893B013BD), // 1.843966568958626e0
    (0x3FFDA9E603DB3285, 0x3E40A3CC98335477), // 1.8539791250833855e0
    (0x3FFDD321F301B460, 0x3E39CEBE40C3061E), // 1.864046048397789e0
    (0x3FFDFC97337B9B5F, 0x3E3DBBC21150A172), // 1.8741676341103e0
    (0x3FFE264614F5A129, 0x3E450E3934FDE7E6), // 1.8843441790323345e0
    (0x3FFE502EE78B3FF6, 0x3E4FDADBBE07FE71), // 1.8945759815869656e0
    (0x3FFE7A51FBC74C83, 0x3E58BBCA5AD6835A), // 1.9048633418176741e0
    (0x3FFEA4AFA2A490DA, 0x3E361428D9F7E139), // 1.9152065613971474e0
    (0x3FFECF482D8E67F1, 0x3E5C28F6CD2CB8AE), // 1.925605943636125e0
    (0x3FFEFA1BEE615A27, 0x3E5DB5DB53E94484), // 1.9360617934922943e0
    (0x3FFF252B376BBA97, 0x3E4E7F799532031A), // 1.9465744175792332e0
    (0x3FFF50765B6E4540, 0x3E575CB1EAC9C7E7), // 1.9571441241754002e0
    (0x3FFF7BFDAD9CBE14, 0x3E5BABA225E6431C), // 1.9677712232331759e0
    (0x3FFFA7C1819E90D8, 0x3E2A3148305BB0F4), // 1.978456026387951e0
    (0x3FFFD3C22B8F71F1, 0x3E573F0776CFD8FB), // 1.9891988469672663e0
];

impl crate::generic::Exp for f64 {
    #[inline]
    fn exp_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp >= 10 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        } else if x_exp < -54 {
            // |x| < 2^-54, so exp(x) = 1 + x + ... rounds to 1
            return 1.0;
        }

        let red = Reduced::exp(x);
        red.eval_to_f64(exp_poly(red.r))
    }

    #[inline]
    fn exp_m1_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp >= 10 {
            if x.is_sign_negative() {
                return -1.0;
            } else {
                return f64::INFINITY;
            }
        } else if x_exp < -12 {
            if x_exp < -53 {
                // exp(x) - 1 = x + x^2 / 2 + ..., where |x| < 2^-53, so
                // x^2 / 2 + ... < 2^-54 * |x|, which is at most half the
                // spacing of `f64` around x, and the result rounds to x
                return x;
            }
            // |x| < 2^-12
            // exp(x) - 1 = x + x^2 * (1/2 + x * Q(x)), where the second term
            // has a relative error of about 3 * 2^-53, so its error relative
            // to the result is at most 2^-64.4
            return x + (x * x) * (0.5 + x * exp_m1_poly(x));
        }

        let red = Reduced::exp(x);
        let k = red.m >> EXP2_TBL_BITS;
        if k < -54 {
            // exp(x) < 2^-54, so exp(x) - 1 rounds to -1
            return -1.0;
        }

        let q = exp_m1_poly(red.r);
        if k > MAX_SCALE_K {
            // exp(x) - 1 rounds to the same value as exp(x), which overflows
            // or is very close to overflowing
            red.eval_to_f64(0.5 + red.r * q)
        } else {
            // exp(x) - 1 = 2^(m / N) * exp(r) - 1 ~= (hi - 1) + lo, where
            // `hi - 1` is calculated exactly as a sum of two `f64` (it is a
            // single `f64` when `-1 <= k <= 52`)
            let (hi, lo) = red.eval_precise(q);
            let d = F64x2::sub11(hi, 1.0);
            d.hi() + (d.lo() + lo)
        }
    }

    #[inline]
    fn exp2_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp >= 11 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        } else if x_exp < -54 {
            // |x * ln(2)| < 2^-54, so exp2(x) rounds to 1
            return 1.0;
        }

        let red = Reduced::exp2(x);
        red.eval_to_f64(exp_poly(red.r))
    }

    #[inline]
    fn exp10_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp >= 9 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        } else if x_exp < -56 {
            // |x * ln(10)| < 2^-54, so exp10(x) rounds to 1
            return 1.0;
        }

        let red = Reduced::exp10(x);
        red.eval_to_f64(exp_poly(red.r))
    }
}

/// Range of `k` for which `2^(m / N) * exp(r)` is calculated with `S`
/// scaled by `2^k`: `S` must not overflow and the small terms must not be
/// subnormal, which would lose precision and can be very slow on some CPUs.
const MIN_SCALE_K: i32 = -960;
const MAX_SCALE_K: i32 = 1023;

/// Returns `P(r)` such that `exp(r) - 1 - r ~= r^2 * P(r)`, with an absolute
/// error of about 2^-64.6, for `|r| <= 0.002711`.
#[inline]
fn exp_poly(r: f64) -> f64 {
    // GENERATE: exp_m1_poly f64 4 -0.002711 0.002711
    const K2: f64 = f64::from_bits(0x3FDFFFFFFFFFFE5A); // 4.999999999999766e-1
    const K3: f64 = f64::from_bits(0x3FC5555555555459); // 1.6666666666665966e-1
    const K4: f64 = f64::from_bits(0x3FA55555C2EEE615); // 4.1666679425770216e-2
    const K5: f64 = f64::from_bits(0x3F81111163F6E1CE); // 8.333335745974415e-3

    let r2 = r * r;
    (K2 + r * K3) + r2 * (K4 + r * K5)
}

/// Returns `Q(r)` such that `exp(r) - 1 - r - r^2 / 2 ~= r^3 * Q(r)`, with an
/// absolute error of about 2^-76.7, which is a relative error of about 2^-68.2
/// in `exp(r) - 1`, for `|r| <= 0.002711`.
///
/// It is more accurate than `exp_poly` because `exp(x) - 1` can be much
/// smaller than `2^(m / N)`.
#[inline]
fn exp_m1_poly(r: f64) -> f64 {
    // GENERATE: exp_m1_poly f64 4 -0.002711 0.002711 2
    const K3: f64 = f64::from_bits(0x3FC55555555554AF); // 1.6666666666666205e-1
    const K4: f64 = f64::from_bits(0x3FA55555555554C8); // 4.1666666666665686e-2
    const K5: f64 = f64::from_bits(0x3F81111156D339E7); // 8.333335363571016e-3
    const K6: f64 = f64::from_bits(0x3F56C16C6D8D4F11); // 1.3888892046524696e-3

    let r2 = r * r;
    (K3 + r * K4) + r2 * (K5 + r * K6)
}

/// Reduced argument: `x = m * ln(2) / N + r`, or the equivalent for other
/// bases.
#[derive(Copy, Clone)]
struct Reduced {
    m: i32,
    /// Multiple of 2^-27 close to `r`, with at most 20 significant bits
    a: f64,
    /// `r - a`
    b: f64,
    /// `r` rounded to `f64`
    r: f64,
}

impl Reduced {
    /// Builds the reduced argument from `m`, `r ~= u + v` (where `|v|` is at
    /// most about 2^-24) and `r` rounded to `f64`.
    #[inline]
    fn new(m: i32, u: f64, v: f64, r: f64) -> Self {
        let a = round_27(u);
        // `u - a` is exact: `a` is a multiple of the ULP of `u` (which is
        // less than 2^-27) and `|u - a| <= |u|`.
        let b = (u - a) + v;
        Self { m, a, b, r }
    }

    /// Reduces `x = m * ln(2) / N + r`.
    ///
    /// `|x|` must be less than 1024.
    #[inline]
    fn exp(x: f64) -> Self {
        let x = x.purify();
        let (mf, m) = round_i32(x * (LOG2_E * TBL_N));
        // |m| < 2^18, so `mf * LN_2_N0` is exact, and so is the subtraction
        // (by Sterbenz lemma, or because `mf = 0`).
        let d = x - mf * LN_2_N0;
        let p = mf * LN_2_N1;
        Self::new(m, d, -p, d - p)
    }

    /// Reduces `x = m / N + r / ln(2)`.
    ///
    /// `|x|` must be less than 2048.
    #[inline]
    fn exp2(x: f64) -> Self {
        let z = x.purify() * TBL_N;
        let (mf, m) = round_i32(z);
        // exact
        let f = z - mf;

        // r = f * ln(2) / N = u + v
        // `f = fh + fl`, where `fh` has 26 bits, so `u` is exact
        // v = fl * X0 + f * X1 = (fl + f * (X1 / X0)) * X0
        let fh = f.split_hi();
        let fl = f - fh;
        let u = fh * (LN_2_X0 / TBL_N);
        let v = (fl + f * (LN_2_X1 / LN_2_X0)) * (LN_2_X0 / TBL_N);
        Self::new(m, u, v, f * (LN_2 / TBL_N))
    }

    /// Reduces `x = m * log10(2) / N + r / ln(10)`.
    ///
    /// `|x|` must be less than 512.
    #[inline]
    fn exp10(x: f64) -> Self {
        let x = x.purify();
        let (mf, m) = round_i32(x * (LOG2_10 * TBL_N));
        // |m| < 2^18, so `mf * LOG10_2_N0` is exact, and so is the
        // subtraction (by Sterbenz lemma, or because `mf = 0`).
        let d = x - mf * LOG10_2_N0;
        let p = mf * LOG10_2_N1;

        // r = (d - p) * ln(10) = u + v
        // `d = dh + dl`, where `dh` has 26 bits, so `u` is exact
        // v = (dl - p) * X0 + (d - p) * X1 = ((dl - p) + (d - p) * (X1 / X0)) * X0
        let dh = d.split_hi();
        let dl = d - dh;
        let u = dh * LN_10_X0;
        let v = ((dl - p) + (d - p) * (LN_10_X1 / LN_10_X0)) * LN_10_X0;
        Self::new(m, u, v, (d - p) * LN_10)
    }

    /// Returns `(hi, lo)` such that `hi + lo ~= 2^(m / N) * exp(r)`, where
    /// `m` replaces `self.m` (only `m mod N` must be the same), with a
    /// relative error of about 2^-64.5 when `exp_poly(r) = p`.
    ///
    /// `hi` is exact and `|lo| < 2^-17 * |hi|`.
    ///
    /// `k = m >> EXP2_TBL_BITS` must be in `[MIN_SCALE_K, MAX_SCALE_K]`.
    #[inline]
    fn eval(&self, m: i32, p: f64) -> (f64, f64) {
        let (s, sh, u) = split_s(m);
        // exact: `sh` has 25 bits and `1 + a` has 28 bits
        let hi = sh * (1.0 + self.a);
        let lo = (sh * self.b + (u + u * self.r)) + (s * (self.r * self.r)) * p;
        (hi, lo)
    }

    /// Like `eval`, but with an absolute error of about 2^-75 relative to
    /// `2^(m / N)` (with `m = self.m`), when `exp_m1_poly(r) = q`.
    ///
    /// `k = m >> EXP2_TBL_BITS` must be in `[MIN_SCALE_K, MAX_SCALE_K]`.
    #[inline]
    fn eval_precise(&self, q: f64) -> (f64, f64) {
        let (s, sh, u) = split_s(self.m);

        // r^2 / 2 = a^2 / 2 + b * (a + b / 2) = w + c, where `a^2 / 2` is exact
        // and `w` is it rounded to a multiple of 2^-27
        let half_a2 = (0.5 * self.a) * self.a;
        let w = round_27(half_a2);
        // `half_a2 - w` is exact
        let c = (half_a2 - w) + self.b * (self.a + 0.5 * self.b);

        // exp(r) = E + (b + c) + r^3 * Q(r), where `E = 1 + a + w` has at
        // most 28 bits (like `1 + a` in `eval`), so `Sh * E` is exact
        // (Sh + U) * exp(r) ~= Sh * E + U * E + S * ((b + c) + r^3 * Q(r))
        let aw = self.a + w;
        let hi = sh * (1.0 + aw);
        let r3 = (self.r * self.r) * self.r;
        let lo = ((u + u * aw) + s * (self.b + c)) + (s * r3) * q;
        (hi, lo)
    }

    /// Returns `2^(m / N) * exp(r)` rounded to `f64`, where `p` is
    /// `exp_poly(r)` or a more accurate approximation.
    #[inline]
    fn eval_to_f64(&self, p: f64) -> f64 {
        let k = self.m >> EXP2_TBL_BITS;
        if (MIN_SCALE_K..=MAX_SCALE_K).contains(&k) {
            let (hi, lo) = self.eval(self.m, p);
            hi + lo
        } else {
            // Calculate `2^(j / N) * exp(r)`, in `[0.99, 2)`, and scale it
            let (hi, lo) = self.eval(self.m & TBL_MASK, p);
            if k > -1022 {
                // The result is normal or overflows, so the scaling is exact
                // (or overflows to infinity).
                scalbn_medium(hi + lo, k)
            } else {
                // The result is subnormal or zero, scale it rounding only
                // once.
                F64x2::fast_add11(hi, lo).scalbn_to_f64(k)
            }
        }
    }
}

/// Returns `(S, Sh, U)`, where `S = 2^k * T`, `Sh` is `S` truncated to 25
/// significant bits, and `U = S * υ`, so `2^(m / N) ~= Sh + U`.
///
/// `k = m >> EXP2_TBL_BITS` must be in `[MIN_SCALE_K, MAX_SCALE_K]`.
#[inline]
fn split_s(m: i32) -> (f64, f64, f64) {
    let (s, upsilon) = exp2_tbl(m);
    // The bits of `T` are truncated like in `EXP2_TBL`.
    let sh = f64::from_bits(s.to_bits() & (u64::MAX << 28));
    (s, sh, s * upsilon)
}

/// Rounds `x` to a multiple of 2^-27 (ties to even).
///
/// `|x|` must be less than 2^24.
#[inline]
fn round_27(x: f64) -> f64 {
    // 1.5 * 2^25, adding it rounds to a multiple of 2^-27.
    const SHIFT: f64 = (3u64 << 24) as f64;

    ((x + SHIFT).purify() - SHIFT).purify()
}

/// Returns `(2^k * T, υ)`, where `m = k * N + j`, `0 <= j < N` and
/// `EXP2_TBL[j] = (T, υ)`. `2^k * T` is `2^(m / N)` rounded to `f64` (see
/// `split_s` for `υ`).
///
/// `2^k * T` must be normal.
#[inline]
pub(crate) fn exp2_tbl(m: i32) -> (f64, f64) {
    let (t_bits, upsilon_bits) = EXP2_TBL[(m & TBL_MASK) as usize];
    let k = m >> EXP2_TBL_BITS;
    debug_assert!((-1022..=1023).contains(&k));
    // Add `k` to the exponent of `T`
    let s = f64::from_bits(t_bits.wrapping_add((k as u64) << 52));
    (s, f64::from_bits(upsilon_bits))
}
