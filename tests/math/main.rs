#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

macro_rules! assert_result_eq {
    ($lhs:expr, $rhs:expr) => {
        let lhs = $lhs;
        let rhs = $rhs;
        if !$crate::ResultEq::result_eq(&lhs, &rhs) {
            panic!(
                "assertion `left == right` (bitwise, with all NaNs equal) failed\n  left: {lhs:?}\n right: {rhs:?}",
            );
        }
    };
}

mod f32;
mod f64;
mod utils;

trait ResultEq {
    fn result_eq(&self, other: &Self) -> bool;
}

impl<T1: ResultEq, T2: ResultEq> ResultEq for (T1, T2) {
    fn result_eq(&self, other: &Self) -> bool {
        self.0.result_eq(&other.0) && self.1.result_eq(&other.1)
    }
}
