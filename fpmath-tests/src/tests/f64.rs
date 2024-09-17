use crate::data::{consume_data, f64 as f64_data, f64::mkfloat};

fn select_threshold(actual: f64, normal_th: f64, subnormal_th: f64) -> f64 {
    if actual == 0.0 || actual.is_subnormal() {
        subnormal_th
    } else {
        normal_th
    }
}

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
            if asin_err > asin_threshold {
                panic!("asin({x:e}) = {actual_asin:e} (error = {asin_err} ULP)");
            }

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            if acos_err > acos_threshold {
                panic!("acos({x:e}) = {actual_acos:e} (error = {acos_err} ULP)");
            }
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
            if asin_err > asin_threshold {
                panic!("asind({x:e}) = {actual_asin:e} (error = {asin_err} ULP)");
            }

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            if acos_err > acos_threshold {
                panic!("acosd({x:e}) = {actual_acos:e} (error = {acos_err} ULP)");
            }
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
            if asin_err > asin_threshold {
                panic!("asinpi({x:e}) = {actual_asin:e} (error = {asin_err} ULP)");
            }

            let acos_threshold = select_threshold(actual_acos, 0.9, 1.9);
            if acos_err > acos_threshold {
                panic!("acospi({x:e}) = {actual_acos:e} (error = {acos_err} ULP)");
            }
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

        if err >= 0.9 {
            panic!("atan({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err >= 0.9 {
            panic!("atand({x:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err >= threshold {
            panic!("atanpi({x:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err > threshold {
            panic!("atan2({x:e}, {y:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err > threshold {
            panic!("atan2d({x:e}, {y:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err > threshold {
            panic!("atan2pi({x:e}, {y:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max atan2pi error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_asinh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_asinh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::asinh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err >= 0.9 {
            panic!("asinh({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err >= 0.9 {
            panic!("acosh({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err >= 0.9 {
            panic!("atanh({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max atanh error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_cbrt() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_cbrt", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::cbrt(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err >= 0.9 {
            panic!("cbrt({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max cbrt error = {max_error}");
    assert!(max_error > 0.5);
}

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
            if exp_err >= exp_threshold {
                panic!("exp({x:e}) = {actual_exp:e} (error = {exp_err} ULP)");
            }

            let expm1_threshold = select_threshold(actual_expm1, 0.9, 1.9);
            if expm1_err >= expm1_threshold {
                panic!("expm1({x:e}) = {actual_expm1:e} (error = {expm1_err} ULP)");
            }
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
        if err > threshold {
            panic!("exp2({x:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err > threshold {
            panic!("exp10({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max exp10 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_hypot() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_hypot", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::hypot(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        if err > threshold {
            panic!("hypot({x:e}, {y:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max hypot error = {max_error}");
    assert!(max_error > 0.49);
}

#[test]
fn test_log() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_log", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::log(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err > 0.9 {
            panic!("log({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err > 0.9 {
            panic!("log_1p({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err > 0.9 {
            panic!("log2({x:e}) = {actual:e} (error = {err} ULP)");
        }
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

        if err > 0.9 {
            panic!("log10({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max log10 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_pow() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_pow", |f64_data::TwoArgData { x, y, expected }| {
        let actual = fpmath::pow(x, y);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        if err > threshold {
            panic!("pow({x:e}, {y:e}) = {actual:e} (error = {err} ULP)");
        }
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
        if err > threshold {
            panic!("powi({x:e}, {y}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max pow error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_sin_cos() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    consume_data(
        "f64_sin_cos",
        |f64_data::SinCosData {
             x,
             expected_sin,
             expected_cos,
         }| {
            let actual_sin1 = fpmath::sin(x);
            let actual_cos1 = fpmath::cos(x);
            let (actual_sin2, actual_cos2) = fpmath::sin_cos(x);

            let sin1_err = expected_sin.calc_error(actual_sin1);
            let sin2_err = expected_sin.calc_error(actual_sin2);
            let cos1_err = expected_cos.calc_error(actual_cos1);
            let cos2_err = expected_cos.calc_error(actual_cos2);

            max_sin1_error = max_sin1_error.max(sin1_err);
            max_sin2_error = max_sin2_error.max(sin2_err);
            max_cos1_error = max_cos1_error.max(cos1_err);
            max_cos2_error = max_cos2_error.max(cos2_err);

            if sin1_err > 0.9 {
                panic!("sin({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)");
            }
            if sin2_err > 0.9 {
                panic!("sin({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)");
            }

            if cos1_err > 0.9 {
                panic!("cos({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)");
            }
            if cos2_err > 0.9 {
                panic!("cos({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)");
            }
        },
    );
    eprintln!("max sin1 error = {max_sin1_error}");
    eprintln!("max sin2 error = {max_sin2_error}");
    eprintln!("max cos1 error = {max_cos1_error}");
    eprintln!("max cos2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_sind_cosd() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    consume_data(
        "f64_sind_cosd",
        |f64_data::SinCosData {
             x,
             expected_sin,
             expected_cos,
         }| {
            let actual_sin1 = fpmath::sind(x);
            let actual_cos1 = fpmath::cosd(x);
            let (actual_sin2, actual_cos2) = fpmath::sind_cosd(x);

            let sin1_err = expected_sin.calc_error(actual_sin1);
            let sin2_err = expected_sin.calc_error(actual_sin2);
            let cos1_err = expected_cos.calc_error(actual_cos1);
            let cos2_err = expected_cos.calc_error(actual_cos2);

            max_sin1_error = max_sin1_error.max(sin1_err);
            max_sin2_error = max_sin2_error.max(sin2_err);
            max_cos1_error = max_cos1_error.max(cos1_err);
            max_cos2_error = max_cos2_error.max(cos2_err);

            let sin1_threshold = select_threshold(actual_sin1, 0.9, 1.9);
            if sin1_err > sin1_threshold {
                panic!("sind({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)");
            }
            let sin2_threshold = select_threshold(actual_sin2, 0.9, 1.9);
            if sin2_err > sin2_threshold {
                panic!("sind({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)");
            }

            let cos1_threshold = select_threshold(actual_cos1, 0.9, 1.9);
            if cos1_err > cos1_threshold {
                panic!("cosd({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)");
            }
            let cos2_threshold = select_threshold(actual_cos2, 0.9, 1.9);
            if cos2_err > cos2_threshold {
                panic!("cosd({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)");
            }
        },
    );
    eprintln!("max sind1 error = {max_sin1_error}");
    eprintln!("max sind2 error = {max_sin2_error}");
    eprintln!("max cosd1 error = {max_cos1_error}");
    eprintln!("max cosd2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

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

            let sin1_err = expected_sin.calc_error(actual_sin1);
            let sin2_err = expected_sin.calc_error(actual_sin2);
            let cos1_err = expected_cos.calc_error(actual_cos1);
            let cos2_err = expected_cos.calc_error(actual_cos2);

            max_sin1_error = max_sin1_error.max(sin1_err);
            max_sin2_error = max_sin2_error.max(sin2_err);
            max_cos1_error = max_cos1_error.max(cos1_err);
            max_cos2_error = max_cos2_error.max(cos2_err);

            let sin1_threshold = select_threshold(actual_sin1, 0.9, 1.9);
            if sin1_err > sin1_threshold {
                panic!("sinpi({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)");
            }
            let sin2_threshold = select_threshold(actual_sin2, 0.9, 1.9);
            if sin2_err > sin2_threshold {
                panic!("sinpi({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)");
            }

            let cos1_threshold = select_threshold(actual_cos1, 0.9, 1.9);
            if cos1_err > cos1_threshold {
                panic!("cospi({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)");
            }
            let cos2_threshold = select_threshold(actual_cos2, 0.9, 1.9);
            if cos2_err > cos2_threshold {
                panic!("cospi({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)");
            }
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

            if sin1_err > 0.9 {
                panic!("sinh({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)");
            }
            if sin2_err > 0.9 {
                panic!("sinh({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)");
            }

            if cos1_err > 0.9 {
                panic!("cosh({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)");
            }
            if cos2_err > 0.9 {
                panic!("cosh({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)");
            }
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
fn test_sqrt() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_sqrt", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::sqrt(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err > 0.5 {
            panic!("sqrt({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max sqrt error = {max_error}");
    assert!(max_error == 0.5);
}

#[test]
fn test_tan() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tan", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tan(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err > 0.9 {
            panic!("tan({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max tan error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tand() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tand", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tand(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        if err > threshold {
            panic!("tand({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max tand error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tanpi() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tanpi", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tanpi(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        if err > threshold {
            panic!("tanpi({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max tanpi error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_tanh() {
    let mut max_error: f64 = 0.0;
    consume_data("f64_tanh", |f64_data::OneArgData { x, expected }| {
        let actual = fpmath::tanh(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        if err > 0.9 {
            panic!("tanh({x:e}) = {actual:e} (error = {err} ULP)");
        }
    });
    eprintln!("max tanh error = {max_error}");
    assert!(max_error > 0.5);
}

pub(super) fn test_round_with(test_f: fn(f64)) {
    use rand::Rng as _;

    let mut rng = crate::data::create_prng();

    for e in -1022..=1023 {
        test_f(mkfloat(0, e, false));
        test_f(mkfloat(0, e, true));
        test_f(mkfloat(u64::MAX, e, false));
        test_f(mkfloat(u64::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.gen::<u64>();
            test_f(mkfloat(m, e, true));
            test_f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=100_000 {
        let arg = arg as f64;
        test_f(arg);
        test_f(-arg);
        test_f(arg + 0.25);
        test_f(-arg + 0.25);
        test_f(arg + 0.5);
        test_f(-arg + 0.5);
        test_f(arg + 0.75);
        test_f(-arg + 0.75);
    }

    for e in 0..=52 {
        for delta in -1000..=1000 {
            let arg = mkfloat(0, e, false) + delta as f64;
            test_f(arg);
            test_f(-arg);
            test_f(arg + 0.25);
            test_f(-arg + 0.25);
            test_f(arg + 0.5);
            test_f(-arg + 0.5);
            test_f(arg + 0.75);
            test_f(-arg + 0.75);
        }
    }
}

#[test]
fn test_round() {
    test_round_with(|arg| {
        let mut expected = dev_mpfr::Mpfr::new();
        expected.set_prec(128);
        expected.set_f64(arg, dev_mpfr::Rnd::N);
        expected.round(None);

        let actual = fpmath::round(arg);

        assert!(
            expected.cmp_f64(actual).is_eq(),
            "round({arg:e}) = {actual:e}",
        );
    });
}

#[test]
fn test_floor() {
    test_round_with(|arg| {
        let mut expected = dev_mpfr::Mpfr::new();
        expected.set_prec(128);
        expected.set_f64(arg, dev_mpfr::Rnd::N);
        expected.floor(None);

        let actual = fpmath::floor(arg);

        assert!(
            expected.cmp_f64(actual).is_eq(),
            "floor({arg:e}) = {actual:e}",
        );
    });
}

#[test]
fn test_ceil() {
    test_round_with(|arg| {
        let mut expected = dev_mpfr::Mpfr::new();
        expected.set_prec(128);
        expected.set_f64(arg, dev_mpfr::Rnd::N);
        expected.ceil(None);

        let actual = fpmath::ceil(arg);

        assert!(
            expected.cmp_f64(actual).is_eq(),
            "ceil({arg:e}) = {actual:e}",
        );
    });
}

#[test]
fn test_trunc() {
    test_round_with(|arg| {
        let mut expected = dev_mpfr::Mpfr::new();
        expected.set_prec(128);
        expected.set_f64(arg, dev_mpfr::Rnd::N);
        expected.trunc(None);

        let actual = fpmath::trunc(arg);

        assert!(
            expected.cmp_f64(actual).is_eq(),
            "trunc({arg:e}) = {actual:e}",
        );
    });
}
