use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_log() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_log", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::log(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log_1p() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_log_1p", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::log_1p(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log_1p({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log_1p error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log2() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_log2", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::log2(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log2({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log10() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_log10", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::log10(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log10({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log10 error = {max_error}");
    assert!(max_error > 0.5);
}
