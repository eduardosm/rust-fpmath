use crate::double::Double;
use crate::scalbn;
use crate::traits::Float as _;

/// Returns `(k, lo, ln_hi)` such as `ln(x) = k * ln(2) + ln_hi + ln(lo)`
#[inline]
pub(crate) fn log_core(x: Double<f32>) -> (f32, Double<f32>, Double<f32>) {
    // GENERATE: ln_table Double<f32> 4
    // LN_TBL[i] = ln(1 + i / 16)       if i < 8
    //           = ln((1 + i / 16) / 2) if i >= 8
    static LN_TBL: [Double<f32>; 16] = [
        Double::new(f32::from_bits(0x00000000), f32::from_bits(0x00000000)), // 0
        Double::new(f32::from_bits(0x3D785186), f32::from_bits(0x2D0B1533)), // 6.0624621816435e-2
        Double::new(f32::from_bits(0x3DF1383B), f32::from_bits(0x3162AF2E)), // 1.1778303565638e-1
        Double::new(f32::from_bits(0x3E2FF983), f32::from_bits(0x32053C9F)), // 1.7185025692666e-1
        Double::new(f32::from_bits(0x3E647FBE), f32::from_bits(0x31735344)), // 2.2314355131421e-1
        Double::new(f32::from_bits(0x3E8B3AE5), f32::from_bits(0x323ABA61)), // 2.7193371548364e-1
        Double::new(f32::from_bits(0x3EA30C5E), f32::from_bits(0x310717B1)), // 3.1845373111853e-1
        Double::new(f32::from_bits(0x3EB9CEBF), f32::from_bits(0x32B5DE80)), // 3.6290549368937e-1
        Double::new(f32::from_bits(0xBE934B10), f32::from_bits(0xB289A6DD)), // -2.8768207245178e-1
        Double::new(f32::from_bits(0xBE7CC8E3), f32::from_bits(0xB1CB3B38)), // -2.4686007793153e-1
        Double::new(f32::from_bits(0xBE549F69), f32::from_bits(0xB26456CF)), // -2.0763936477824e-1
        Double::new(f32::from_bits(0xBE2DFA03), f32::from_bits(0xB1B543DB)), // -1.6989903679540e-1
        Double::new(f32::from_bits(0xBE08BC74), f32::from_bits(0xB089F91F)), // -1.3353139262452e-1
        Double::new(f32::from_bits(0xBDC99AF2), f32::from_bits(0xB1EACA4C)), // -9.8440072813253e-2
        Double::new(f32::from_bits(0xBD842CC5), f32::from_bits(0xB1ACF1D0)), // -6.4538521137571e-2
        Double::new(f32::from_bits(0xBD020AEC), f32::from_bits(0xB09E7444)), // -3.1748698314580e-2
    ];

    // GENERATE: ln_lo_scale_table Double<f32> 4
    // LN_LO_SCALE_TBL[i] = 1 / (1 + i / 16)
    static LN_LO_SCALE_TBL: [Double<f32>; 16] = [
        Double::new(f32::from_bits(0x3F800000), f32::from_bits(0x00000000)), // 1.0000000000000e0
        Double::new(f32::from_bits(0x3F70F0F0), f32::from_bits(0x3370F0F1)), // 9.4117647058824e-1
        Double::new(f32::from_bits(0x3F638E38), f32::from_bits(0x33638E39)), // 8.8888888888889e-1
        Double::new(f32::from_bits(0x3F579435), f32::from_bits(0x33650D79)), // 8.4210526315789e-1
        Double::new(f32::from_bits(0x3F4CCCCC), f32::from_bits(0x334CCCCD)), // 8.0000000000000e-1
        Double::new(f32::from_bits(0x3F430C30), f32::from_bits(0x33430C31)), // 7.6190476190476e-1
        Double::new(f32::from_bits(0x3F3A2E8B), f32::from_bits(0x3322E8BA)), // 7.2727272727273e-1
        Double::new(f32::from_bits(0x3F321642), f32::from_bits(0x3348590B)), // 6.9565217391304e-1
        Double::new(f32::from_bits(0x3F2AAAAA), f32::from_bits(0x332AAAAB)), // 6.6666666666667e-1
        Double::new(f32::from_bits(0x3F23D70A), f32::from_bits(0x3275C28F)), // 6.4000000000000e-1
        Double::new(f32::from_bits(0x3F1D89D8), f32::from_bits(0x331D89D9)), // 6.1538461538462e-1
        Double::new(f32::from_bits(0x3F17B425), f32::from_bits(0x336D097B)), // 5.9259259259259e-1
        Double::new(f32::from_bits(0x3F124924), f32::from_bits(0x33124925)), // 5.7142857142857e-1
        Double::new(f32::from_bits(0x3F0D3DCB), f32::from_bits(0x310D3DCB)), // 5.5172413793103e-1
        Double::new(f32::from_bits(0x3F088888), f32::from_bits(0x33088889)), // 5.3333333333333e-1
        Double::new(f32::from_bits(0x3F042108), f32::from_bits(0x32842108)), // 5.1612903225806e-1
    ];

    // x = 2^k * m, with 1 <= m < 2
    let k = x.hi().exponent();
    let m = x.hi().mant();

    // ln(m) = ln(hi) + ln(lo + 1)
    let hi = m >> (24 - 5) & 0xF;
    let hi_f = f32::from_bits(0x3F80_0000 | (hi << (24 - 5)));
    let hi = hi as usize;
    let scale = scalbn(1.0, (-k).into());
    let lo = (x.pmul1(scale) - hi_f) * LN_LO_SCALE_TBL[hi];
    let ln_hi = LN_TBL[hi];

    // The upper half of `LN_TBL` holds `ln((1 + i / 16) / 2)` instead of
    // `ln(1 + i / 16)` to avoid cancellation when `x` is close to 1.
    let k = if hi >= 8 { k + 1 } else { k };

    (f32::from(k), lo, ln_hi)
}
