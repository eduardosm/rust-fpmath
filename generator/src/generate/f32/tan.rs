use super::super::sollya;

pub(in super::super) fn gen_tan_poly() -> String {
    let mut out = String::new();

    let f = "tan(x) / x - 1";
    let poly_i = [2, 4, 6, 8];
    let range0 = 0.393; // ~= π/8
    let range = (-range0, range0);

    sollya::run_and_render_remez_f32(f, range, &poly_i, 1, "K", &mut out);

    out
}
