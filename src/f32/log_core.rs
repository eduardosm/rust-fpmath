use crate::traits::Float as _;

// GENERATE: ln_table f32 f64 4
// LN_TBL[i] = -ln(LN_LO_SCALE_TBL[i])
static LN_TBL: [f64; 17] = [
    f64::from_bits(0x8000000000000000), // -0e0
    f64::from_bits(0x3FAF0A30A01162A7), // 6.062461809114455e-2
    f64::from_bits(0x3FBE27074E2AF2E8), // 1.1778302820580289e-1
    f64::from_bits(0x3FC5FF3060A793D5), // 1.7185024947607866e-1
    f64::from_bits(0x3FCC8FF7A79A9A26), // 2.2314353641304868e-1
    f64::from_bits(0x3FD1675C97ABA611), // 2.7193369685719043e-1
    f64::from_bits(0x3FD4618BA21C5ECA), // 3.1845370131621265e-1
    f64::from_bits(0x3FD739D7E2BBD00A), // 3.6290547506291715e-1
    f64::from_bits(0x3FD9F323CCBF9854), // 4.0546507830584244e-1
    f64::from_bits(0x3FDC8FF7DF9A9A26), // 4.4628712498016154e-1
    f64::from_bits(0x3FDF128F37AF06F9), // 4.855077785287985e-1
    f64::from_bits(0x3FE0BE72E0252A83), // 5.232481363139673e-1
    f64::from_bits(0x3FE1E85F467040D9), // 5.596157432319401e-1
    f64::from_bits(0x3FE307D7354F10BE), // 5.947071114719831e-1
    f64::from_bits(0x3FE41D8FCC4672BB), // 6.286086072683114e-1
    f64::from_bits(0x3FE52A2D365BC5AF), // 6.613985120476878e-1
    f64::from_bits(0x3FE62E42FEFA39EF), // 6.931471805599453e-1
];

// GENERATE: ln_lo_scale_table f32 4
// LN_LO_SCALE_TBL[i] = 1 / (1 + i / 16)
static LN_LO_SCALE_TBL: [f32; 17] = [
    f32::from_bits(0x3F800000), // 1e0
    f32::from_bits(0x3F70F0F1), // 9.411765e-1
    f32::from_bits(0x3F638E39), // 8.888889e-1
    f32::from_bits(0x3F579436), // 8.4210527e-1
    f32::from_bits(0x3F4CCCCD), // 8e-1
    f32::from_bits(0x3F430C31), // 7.619048e-1
    f32::from_bits(0x3F3A2E8C), // 7.2727275e-1
    f32::from_bits(0x3F321643), // 6.956522e-1
    f32::from_bits(0x3F2AAAAB), // 6.666667e-1
    f32::from_bits(0x3F23D70A), // 6.4e-1
    f32::from_bits(0x3F1D89D9), // 6.1538464e-1
    f32::from_bits(0x3F17B426), // 5.925926e-1
    f32::from_bits(0x3F124925), // 5.714286e-1
    f32::from_bits(0x3F0D3DCB), // 5.5172414e-1
    f32::from_bits(0x3F088889), // 5.3333336e-1
    f32::from_bits(0x3F042108), // 5.16129e-1
    f32::from_bits(0x3F000000), // 5e-1
];

#[inline]
pub(crate) fn log_core_f32(x: f32, edelta: i16) -> (f64, f64, f64) {
    // Split x * 2^edelta = 2^k * m
    // m = hi * (1 + lo)

    let m = x.mant();
    let hi = ((m >> (24 - 5 - 1) & 0x1F).div_ceil(2)) as usize;
    let ln_hi = LN_TBL[hi];
    // lo = m / hi - 1
    let lo = f64::from(x.set_exp(0)) * f64::from(LN_LO_SCALE_TBL[hi]) - 1.0;

    let k = x.exponent() + edelta;

    (f64::from(k), lo, ln_hi)
}

#[inline]
pub(crate) fn log_core_f64(x: f64, edelta: i16) -> (f64, f64, f64) {
    // Split x * 2^edelta = 2^k * m
    // m = hi * (1 + lo)

    let m = x.mant();
    let hi = ((m >> (53 - 5 - 1) & 0x1F).div_ceil(2)) as usize;
    let ln_hi = LN_TBL[hi];
    // lo = m / hi - 1
    let lo = x.set_exp(0) * f64::from(LN_LO_SCALE_TBL[hi]) - 1.0;

    let k = x.exponent() + edelta;

    (f64::from(k), lo, ln_hi)
}
