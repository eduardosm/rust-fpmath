use super::select_threshold;
use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_pow() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_pow", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::pow(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "pow({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max pow error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_powi() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_powi", |f64_data::powi::Data { x, y, expected }| {
        let actual = fpmath::powi(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "powi({x:e}, {y}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max pow error = {max_error}");
    assert!(max_error > 0.5);
}
