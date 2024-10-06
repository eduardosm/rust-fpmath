use super::select_threshold;
use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_exp() {
    let mut max_exp_error: f64 = 0.0;
    let mut max_expm1_error: f64 = 0.0;
    consume_data(
        "f64_exp",
        |f64_data::exp::ExpExpM1Data {
             x,
             expected_exp,
             expected_expm1,
         }| {
            let actual_exp = fpmath::exp(x);
            let actual_expm1 = fpmath::exp_m1(x);

            let exp_err = expected_exp.calc_error(actual_exp);
            let expm1_err = expected_expm1.calc_error(actual_expm1);

            max_exp_error = max_exp_error.max(exp_err);
            max_expm1_error = max_expm1_error.max(expm1_err);

            let exp_threshold = select_threshold(actual_exp, 0.9, 1.9);
            assert!(
                exp_err < exp_threshold,
                "exp({x:e}) = {actual_exp:e} (error = {exp_err} ULP)",
            );

            let expm1_threshold = select_threshold(actual_expm1, 0.9, 1.9);
            assert!(
                expm1_err < expm1_threshold,
                "expm1({x:e}) = {actual_expm1:e} (error = {expm1_err} ULP)",
            );
        },
    );
    eprintln!("max exp error = {max_exp_error}");
    eprintln!("max expm1 error = {max_expm1_error}");
    assert!(max_exp_error > 0.5);
    assert!(max_expm1_error > 0.5);
}

#[test]
fn test_exp2() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_exp2", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::exp2(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "exp2({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max exp2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_exp10() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_exp10", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::exp10(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "exp10({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max exp10 error = {max_error}");
    assert!(max_error > 0.5);
}
