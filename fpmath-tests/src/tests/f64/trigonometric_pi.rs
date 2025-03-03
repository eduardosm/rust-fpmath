use super::{purify, purify2, select_threshold};
use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_sinpi_cospi() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    consume_data(
        "f64_sinpi_cospi",
        |f64_data::SinCosData {
             x,
             expected_sin,
             expected_cos,
         }| {
            let actual_sin1 = fpmath::sinpi(x);
            let actual_cos1 = fpmath::cospi(x);
            let (actual_sin2, actual_cos2) = fpmath::sinpi_cospi(x);
            assert_eq!(purify(fpmath::sinpi(-x)), purify(-actual_sin1));
            assert_eq!(purify(fpmath::cospi(-x)), purify(actual_cos1));
            assert_eq!(
                purify2(fpmath::sinpi_cospi(-x)),
                purify2((-actual_sin2, actual_cos2))
            );

            let sin1_err = expected_sin.calc_error(actual_sin1);
            let sin2_err = expected_sin.calc_error(actual_sin2);
            let cos1_err = expected_cos.calc_error(actual_cos1);
            let cos2_err = expected_cos.calc_error(actual_cos2);

            max_sin1_error = max_sin1_error.max(sin1_err);
            max_sin2_error = max_sin2_error.max(sin2_err);
            max_cos1_error = max_cos1_error.max(cos1_err);
            max_cos2_error = max_cos2_error.max(cos2_err);

            let sin1_threshold = select_threshold(actual_sin1, 0.9, 1.9);
            assert!(
                sin1_err < sin1_threshold,
                "sinpi({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
            );
            let sin2_threshold = select_threshold(actual_sin2, 0.9, 1.9);
            assert!(
                sin2_err < sin2_threshold,
                "sinpi({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
            );

            let cos1_threshold = select_threshold(actual_cos1, 0.9, 1.9);
            assert!(
                cos1_err < cos1_threshold,
                "cospi({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
            );
            let cos2_threshold = select_threshold(actual_cos2, 0.9, 1.9);
            assert!(
                cos2_err < cos2_threshold,
                "cospi({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
            );
        },
    );
    eprintln!("max sinpi1 error = {max_sin1_error}");
    eprintln!("max sinpi2 error = {max_sin2_error}");
    eprintln!("max cospi1 error = {max_cos1_error}");
    eprintln!("max cospi2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_tanpi() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tanpi", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tanpi(x);
        assert_eq!(purify(fpmath::tanpi(-x)), purify(-actual));

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "tanpi({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max tanpi error = {max_error}");
    assert!(max_error > 0.5);
}
