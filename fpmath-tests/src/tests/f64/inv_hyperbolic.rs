use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_asinh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_asinh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::asinh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "asinh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max asinh error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_acosh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_acosh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::acosh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "acosh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max acosh error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atanh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atanh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::atanh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atanh error = {max_error}");
    assert!(max_error > 0.5);
}
