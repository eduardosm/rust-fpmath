use crate::data::{consume_data, f64 as f64_data};

#[test]
fn test_tgamma() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tgamma", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tgamma(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = if x < 0.5 { 1.9 } else { 0.9 };
        assert!(
            err < threshold,
            "tgamma({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max tgamma error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_lgamma() {
    let mut max_error: f64 = 0.0;
    consume_data(
        "f64_lgamma",
        |f64_data::gamma::LgammaData {
             x,
             expected,
             expected_sign,
         }| {
            let (actual, actual_sign) = fpmath::lgamma(x);

            let err = expected.calc_error(actual);
            max_error = max_error.max(err);

            let threshold = if (-5.0..=-2.0).contains(&x) {
                // FIXME
                50.0
            } else if (0.5..=7.0).contains(&x) {
                1.5
            } else {
                1.9
            };
            assert_eq!(expected_sign, actual_sign);
            assert!(
                err < threshold,
                "lgamma({x:e}) = {actual:e} (error = {err} ULP)",
            );
        },
    );
    eprintln!("max lgamma error = {max_error}");
    assert!(max_error > 0.5);
}
