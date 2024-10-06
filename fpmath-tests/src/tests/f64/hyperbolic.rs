use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_sinh_cosh() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    consume_data(
        "f64_sinh_cosh",
        |f64_data::SinCosData {
             x,
             expected_sin,
             expected_cos,
         }| {
            let actual_sin1 = fpmath::sinh(x);
            let actual_cos1 = fpmath::cosh(x);
            let (actual_sin2, actual_cos2) = fpmath::sinh_cosh(x);

            let sin1_err = expected_sin.calc_error(actual_sin1);
            let sin2_err = expected_sin.calc_error(actual_sin2);
            let cos1_err = expected_cos.calc_error(actual_cos1);
            let cos2_err = expected_cos.calc_error(actual_cos2);

            max_sin1_error = max_sin1_error.max(sin1_err);
            max_sin2_error = max_sin2_error.max(sin2_err);
            max_cos1_error = max_cos1_error.max(cos1_err);
            max_cos2_error = max_cos2_error.max(cos2_err);

            assert!(
                sin1_err < 0.9,
                "sinh({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
            );
            assert!(
                sin2_err < 0.9,
                "sinh({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
            );

            assert!(
                cos1_err < 0.9,
                "cosh({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
            );
            assert!(
                cos2_err < 0.9,
                "cosh({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
            );
        },
    );
    eprintln!("max sinh1 error = {max_sin1_error}");
    eprintln!("max sinh2 error = {max_sin2_error}");
    eprintln!("max cosh1 error = {max_cos1_error}");
    eprintln!("max cosh2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_tanh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tanh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tanh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "tanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max tanh error = {max_error}");
    assert!(max_error > 0.5);
}
