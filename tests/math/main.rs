#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

mod f32;
mod f64;

fn create_prng() -> impl rand::Rng {
    use rand::SeedableRng as _;
    rand_pcg::Pcg64::seed_from_u64(0x985A_9231_A004_6A3D)
}
