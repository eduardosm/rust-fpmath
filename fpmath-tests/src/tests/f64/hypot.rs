use super::select_threshold;
use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_hypot() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_hypot", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::hypot(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "hypot({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max hypot error = {max_error}");
    assert!(max_error > 0.49);
}
