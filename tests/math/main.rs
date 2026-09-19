#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

macro_rules! assert_total_eq {
    ($lhs:expr, $rhs:expr) => {
        let lhs = $lhs;
        let rhs = $rhs;
        if !$crate::TotalEq::total_eq(&lhs, &rhs) {
            panic!(
                "assertion `left == right` (using totalOrder) failed\n  left: {lhs:?}\n right: {rhs:?}",
            );
        }
    };
}

mod f32;
mod f64;

fn create_prng() -> impl rand::Rng {
    use rand::SeedableRng as _;
    rand_pcg::Pcg64::seed_from_u64(0x985A_9231_A004_6A3D)
}

trait TotalEq {
    fn total_eq(&self, other: &Self) -> bool;
}

impl<T1: TotalEq, T2: TotalEq> TotalEq for (T1, T2) {
    fn total_eq(&self, other: &Self) -> bool {
        self.0.total_eq(&other.0) && self.1.total_eq(&other.1)
    }
}
