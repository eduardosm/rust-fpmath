use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_cbrt() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_cbrt", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::cbrt(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "cbrt({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max cbrt error = {max_error}");
    assert!(max_error > 0.5);
}
