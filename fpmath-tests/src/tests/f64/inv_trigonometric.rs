use super::select_threshold;
use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_asin_acos() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    consume_data(
        "f64_asin_acos",
        |f64_data::asin_acos::Data {
             x,
             expected_asin,
             expected_acos,
         }| {
            let actual_asin = fpmath::asin(x);
            let actual_acos = fpmath::acos(x);

            let asin_err = expected_asin.calc_error(actual_asin);
            let acos_err = expected_acos.calc_error(actual_acos);

            max_asin_error = max_asin_error.max(asin_err);
            max_acos_error = max_acos_error.max(acos_err);

            let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
            assert!(
                asin_err < asin_threshold,
                "asin({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
            );

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            assert!(
                acos_err < acos_threshold,
                "acos({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
            );
        },
    );
    eprintln!("max asin error = {max_asin_error}");
    eprintln!("max acos error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

#[test]
fn test_asind_acosd() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    consume_data(
        "f64_asind_acosd",
        |f64_data::asin_acos::Data {
             x,
             expected_asin,
             expected_acos,
         }| {
            let actual_asin = fpmath::asind(x);
            let actual_acos = fpmath::acosd(x);

            let asin_err = expected_asin.calc_error(actual_asin);
            let acos_err = expected_acos.calc_error(actual_acos);

            max_asin_error = max_asin_error.max(asin_err);
            max_acos_error = max_acos_error.max(acos_err);

            let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
            assert!(
                asin_err < asin_threshold,
                "asind({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
            );

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            assert!(
                acos_err < acos_threshold,
                "acosd({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
            );
        },
    );
    eprintln!("max asind error = {max_asin_error}");
    eprintln!("max acosd error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

#[test]
fn test_asinpi_acospi() {
    let mut max_asin_error: f64 = 0.0;
    let mut max_acos_error: f64 = 0.0;
    consume_data(
        "f64_asinpi_acospi",
        |f64_data::asin_acos::Data {
             x,
             expected_asin,
             expected_acos,
         }| {
            let actual_asin = fpmath::asinpi(x);
            let actual_acos = fpmath::acospi(x);

            let asin_err = expected_asin.calc_error(actual_asin);
            let acos_err = expected_acos.calc_error(actual_acos);

            max_asin_error = max_asin_error.max(asin_err);
            max_acos_error = max_acos_error.max(acos_err);

            let asin_threshold = select_threshold(actual_asin, 0.9, 1.9);
            assert!(
                asin_err < asin_threshold,
                "asinpi({x:e}) = {actual_asin:e} (error = {asin_err} ULP)",
            );

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            assert!(
                acos_err < acos_threshold,
                "acospi({x:e}) = {actual_acos:e} (error = {acos_err} ULP)",
            );
        },
    );
    eprintln!("max asinpi error = {max_asin_error}");
    eprintln!("max acospi error = {max_acos_error}");
    assert!(max_asin_error > 0.5);
    assert!(max_acos_error > 0.5);
}

#[test]
fn test_atan() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atan", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::atan(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atan({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atan error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atand() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atand", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::atand(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atand({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atand error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atanpi() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atanpi", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::atanpi(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atanpi({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atanpi error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atan2", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::atan2(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atan2({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2d() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atan2d", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::atan2d(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "atan2d({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2d error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_atan2pi() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_atan2pi", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::atan2pi(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.95, 1.9);
        assert!(
            err < threshold,
            "atan2pi({x:e}, {y:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max atan2pi error = {max_error}");
    assert!(max_error > 0.5);
}
