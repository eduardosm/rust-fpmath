use crate::sollya;

pub(crate) fn gen_sin_poly() {
    let f = "sin(x) / x - 1";
    let poly_i = [2, 4, 6, 99];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 1, "K");
}

pub(crate) fn gen_cos_poly() {
    let f = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = [4, 6, 8, 99];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 0, "K");
}
