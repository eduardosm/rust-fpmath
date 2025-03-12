// 2 / π ~= sum(FRAC_2_PI_LARGE[i] * 2^(-24 * (i + 1)))
// GENERATE: common::reduce_pi_2_large::FRAC_2_PI_LARGE
const FRAC_2_PI_LARGE: &[u32] = &[
    0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62, 0x95993C, 0x439041, 0xFE5163,
    0xABDEBB, 0xC561B7, 0x246E3A, 0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
    0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41, 0x3991D6, 0x398353, 0x39F49C,
    0x845F8B, 0xBDF928, 0x3B1FF8, 0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
    0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5, 0xF17B3D, 0x0739F7, 0x8A5292,
    0xEA6BFB, 0x5FB11F, 0x8D5D08, 0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
    0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880, 0x4D7327, 0x310606, 0x1556CA,
    0x73A8C9, 0x60E27B, 0xC08C6B,
];

// π / 2 ~= sum(FRAC_PI_2_MEDIUM[i] * 2^(1 - 24 * (i + 1)))
// GENERATE: common::reduce_pi_2_large::FRAC_PI_2_MEDIUM
const FRAC_PI_2_MEDIUM: &[u32] = &[
    0xC90FDA, 0xA22168, 0xC234C4, 0xC6628B, 0x80DC1C, 0xD12902, 0x4E088A, 0x67CC74,
];

// Based on __rem_pio2_large from musl libc, which has been
// translated from C to Rust (trying to be idiomatic). Floating
// point arithmetic has been replaced with integer or fixed
// point arithmetic.

// ORIGINAL COPYRIGHT NOTICE
//
// origin: FreeBSD /usr/src/lib/msun/src/k_rem_pio2.c
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================

/// Common part of π/2 reduction
///
/// Inputs:
/// * `x` is the input splitted in 24-bit chunks
/// * `e0` is the exponent of `x[0]` (`x[0]*2^e0` is the input truncated to 24
///   bits).
/// * `jk + 1` is the initial number of needed terms of `FRAC_2_PI_LARGE`, which
///   is 4 for `f32` and 5 for `f64`.
///
/// Returns `(ih, jz, n, qe)`.
pub(crate) fn reduce_pi_2_large(
    x: &[u32],
    e0: i16,
    jk: usize,
    qp: &mut [u64; 20],
) -> (u32, usize, u8, i16) {
    assert!(e0 >= -24);

    // jp + 1 is the number of needed terms of FRAC_PI_2_MEDIUM,
    // which is equal to jk + 1.
    let jp = jk;

    // jv is the index pointing to the position of FRAC_2_PI_LARGE
    // needed to calculate q.
    let jv = ((e0 - 3) / 24).max(0);
    // Exponent of q[0], such as the exponent of q[i] is q0 - 24 * i.
    // Note that q0 < 3
    let mut q0 = e0 - 24 * (jv + 1);
    debug_assert!(matches!(q0, -48..=2));
    let jv = jv as usize;

    // Calculate q[0..=jk] (array of chunks representing the product
    // of x and 2/π)
    let mut q: [u64; 20] = [0; 20];
    for i in 0..=jk {
        let mut fw = 0;
        for j in 0..x.len().min(jv + i + 1) {
            fw += u64::from(x[j]) * u64::from(FRAC_2_PI_LARGE[jv + i - j]);
        }
        q[i] = fw;
    }

    // jz + 1 is the number of used terms of FRAC_2_PI_LARGE.
    let mut jz = jk;
    // q splitted in 24-bit chunks
    let mut iq: [u32; 20] = [0; 20];
    let mut z: u64;
    let mut n: u32;
    // if ih > 0, q >= 0.5
    let mut ih: u32;

    loop {
        // Distill q[] into iq[] reversingly
        z = q[jz];
        for (iq_i, q_i) in iq[..jz].iter_mut().zip(q[..jz].iter().rev()) {
            *iq_i = (z & 0x00FF_FFFF) as u32;
            z = *q_i + (z >> 24);
        }

        // Calculate n and remove integer part from z
        if q0 >= 0 {
            n = ((z << q0) % 8) as u32;
            z = 0;
        } else {
            // q0 >= -48 && q0 < 0
            n = ((z >> -q0) % 8) as u32;
            z &= !(u64::MAX << (-q0));
        }

        ih = 0;
        if q0 > 0 {
            // Need iq[jz-1] to determine n
            n += iq[jz - 1] >> (24 - q0);
            iq[jz - 1] &= (1 << (24 - q0)) - 1;
            ih = iq[jz - 1] >> (23 - q0);
            debug_assert!(ih <= 1);
        } else if q0 == 0 {
            ih = iq[jz - 1] >> 23;
            debug_assert!(ih <= 1);
        } else if z >= (1 << -(q0 + 1)) {
            // q0 >= -48 && q0 < 0
            // z >= 0.5
            ih = 2;
        }

        if ih > 0 {
            // q > 0.5
            n += 1;
            let mut carry = false;
            for iq_i in iq[0..jz].iter_mut() {
                // Compute 1-q
                let j = *iq_i;
                if !carry {
                    if j != 0 {
                        carry = true;
                        *iq_i = 0x0100_0000 - j;
                    }
                } else {
                    *iq_i = 0x00FF_FFFF - j;
                }
            }
            match q0 {
                // Rare case: chance is 1 in 12
                1 => iq[jz - 1] &= 0x007F_FFFF,
                2 => iq[jz - 1] &= 0x003F_FFFF,
                _ => {}
            }
            if ih == 2 {
                debug_assert!(q0 < 0);
                z = (1 << -q0) - z;
                if carry {
                    z -= 1;
                }
            }
        }

        // Check if recomputation is needed
        if z == 0 {
            let mut j = 0;
            for i in (jk..jz).rev() {
                j |= iq[i];
            }
            if j == 0 {
                // Recomputation needed
                // k is the number of needed terms
                let k = iq[..jk].iter().rev().take_while(|&&iq_i| iq_i == 0).count() + 1;

                for i in (jz + 1)..=(jz + k) {
                    // Add q[jz+1] to q[jz+k]
                    let mut fw = 0;
                    for j in 0..x.len() {
                        fw += u64::from(x[j]) * u64::from(FRAC_2_PI_LARGE[jv + i - j]);
                    }
                    q[i] = fw;
                }
                jz += k;
                continue;
            }
        }

        break;
    }

    // Chop off zero terms
    if z == 0 {
        jz -= 1;
        q0 -= 24;
        while iq[jz] == 0 {
            jz -= 1;
            q0 -= 24;
        }
    } else {
        // Break z into 24-bit if necessary
        if z >= 0x0100_0000 {
            iq[jz] = (z & 0x00FF_FFFF) as u32;
            jz += 1;
            q0 += 24;
            iq[jz] = (z >> 24) as u32;
        } else {
            iq[jz] = z as u32;
        }
    }

    // Compute FRAC_PI_2_MEDIUM[0..=jp] * q[jz..=0]
    let qe = q0 - 23;
    for i in (0..=jz).rev() {
        let mut acc = 0;
        for k in 0..=jp.min(jz - i) {
            acc += u64::from(FRAC_PI_2_MEDIUM[k]) * u64::from(iq[i + k]);
        }
        qp[jz - i] = acc;
    }

    (ih, jz, n as u8, qe)
}
