mod cbrt;
mod exp;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod round;
mod sqrt;
mod trigonometric_deg;
mod trigonometric_pi;
mod trigonometric_rad;

fn select_threshold(actual: f64, normal_th: f64, subnormal_th: f64) -> f64 {
    if actual == 0.0 || actual.is_subnormal() {
        subnormal_th
    } else {
        normal_th
    }
}
