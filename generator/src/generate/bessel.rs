use std::fmt::Write as _;

use rug::ops::Pow as _;

use super::{
    CoeffKinds, FloatKind, arg_utils, julia, render_const, render_const_dec_value,
    render_const_value,
};

/// Precision of the auxiliary calculations, which is much higher than the one
/// of the results because some of them have large cancellations (see
/// `Func::anchored_taylor`).
const PREC: u32 = 2048;

/// Number of decimal digits of the values passed to Julia.
const JULIA_DIGITS: usize = 150;

/// Number of terms of the Taylor series that represent the functions in the
/// intervals of the tables (they are much more accurate than the fitted
/// polynomials).
const TAYLOR_TERMS: usize = 80;

/// A Bessel function of the first (`J`) or second (`Y`) kind, of order 0 or 1.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Func {
    J0,
    J1,
    Y0,
    Y1,
}

impl std::str::FromStr for Func {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "j0" => Ok(Self::J0),
            "j1" => Ok(Self::J1),
            "y0" => Ok(Self::Y0),
            "y1" => Ok(Self::Y1),
            _ => Err("invalid Bessel function"),
        }
    }
}

impl Func {
    fn name(self) -> &'static str {
        match self {
            Self::J0 => "j0",
            Self::J1 => "j1",
            Self::Y0 => "y0",
            Self::Y1 => "y1",
        }
    }

    fn order(self) -> u32 {
        match self {
            Self::J0 | Self::Y0 => 0,
            Self::J1 | Self::Y1 => 1,
        }
    }

    /// Returns the function of the first kind of the same order.
    fn first_kind(self) -> Self {
        match self {
            Self::J0 | Self::Y0 => Self::J0,
            Self::J1 | Self::Y1 => Self::J1,
        }
    }

    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = x.clone();
        match self {
            Self::J0 => x.j0(),
            Self::J1 => x.j1(),
            Self::Y0 => x.y0(),
            Self::Y1 => x.y1(),
        }
    }

    /// Calculates the derivative of the function, with `C0'(x) = -C1(x)` and
    /// `C1'(x) = C0(x) - C1(x) / x`.
    fn eval_deriv(self, x: &rug::Float) -> rug::Float {
        match self {
            Self::J0 => -x.clone().j1(),
            Self::J1 => x.clone().j0() - x.clone().j1() / x,
            Self::Y0 => -x.clone().y1(),
            Self::Y1 => x.clone().y0() - x.clone().y1() / x,
        }
    }

    /// The `k`-th positive zero of the function (counting from 1) is close to
    /// `(k + zero_offset) * π`.
    fn zero_offset(self) -> f64 {
        match self {
            Self::J0 => -0.25,
            Self::J1 => 0.25,
            Self::Y0 => -0.75,
            Self::Y1 => -0.25,
        }
    }

    /// Calculates the `k`-th positive zero of the function (counting from 1).
    fn zero(self, k: u32) -> Result<rug::Float, String> {
        let pi = rug::Float::with_val(PREC, rug::float::Constant::Pi);
        // McMahon's expansion, with b = (k + zero_offset) * π and μ = 4 * ν^2:
        // b - (μ - 1) / (8 * b) - ...
        let b = pi * (f64::from(k) + self.zero_offset());
        let mu = 4 * self.order().pow(2);
        let approx: rug::Float = b.clone() - (f64::from(mu) - 1.0) / (8 * b);

        // Refine the zero with Newton's method
        let mut x = approx.clone();
        for _ in 0..20 {
            let step = self.eval(&x) / self.eval_deriv(&x);
            x -= step;
        }
        let residual = self.eval(&x);
        if !residual.is_zero() && residual.get_exp().unwrap() > -(PREC as i32 - 64) {
            return Err(format!("zero {k} of {} did not converge", self.name()));
        }
        if (x.clone() - &approx).abs() > 0.1 {
            return Err(format!(
                "zero {k} of {} converged to {x:e}, far from {approx:e}",
                self.name()
            ));
        }
        Ok(x)
    }

    /// Returns the first `n` coefficients of the Taylor series of the function
    /// around `c`, in powers of `d = x - c`.
    ///
    /// `y = sum(a[k] * d^k)` is a solution of
    /// `x^2 * y'' + x * y' + (x^2 - ν^2) * y = 0`, so
    ///
    /// `c^2 * (k + 2) * (k + 1) * a[k+2] = -(p * a[k+1] + q * a[k] + 2 * c * a[k-1] + a[k-2])`
    ///
    /// with `p = c * (k + 1) * (2 * k + 1)` and `q = k^2 + c^2 - ν^2`.
    ///
    /// The errors of this recurrence grow like the coefficients of the
    /// functions of the second kind (`c^-k`), which are much larger than the
    /// ones of the functions of the first kind, but their effect on the sum
    /// is negligible with the precision used here when `|d| < c`.
    fn taylor(self, c: &rug::Float, n: usize) -> Vec<rug::Float> {
        let nu2 = self.order().pow(2);
        let c2 = rug::Float::with_val(PREC, c.square_ref());
        let mut a = vec![self.eval(c), self.eval_deriv(c)];
        for k in 0..n.saturating_sub(2) {
            let ku = k as u32;
            let mut s = rug::Float::with_val(PREC, &a[k + 1] * c) * ((ku + 1) * (2 * ku + 1));
            s += (c2.clone() + ku * ku - nu2) * &a[k];
            if k >= 1 {
                s += rug::Float::with_val(PREC, &a[k - 1] * c) * 2u32;
            }
            if k >= 2 {
                s += &a[k - 2];
            }
            let den = c2.clone() * ((ku + 2) * (ku + 1));
            a.push(-s / den);
        }
        a.truncate(n);
        a
    }

    /// Returns the first `n` coefficients of the Taylor series around `c`, in
    /// powers of `d = x - c`, of the function without its singularity at
    /// zero, which makes it entire:
    ///
    /// * `Y0(x) - 2/π * J0(x) * ln(x / c)`
    /// * `x * Y1(x) - 2/π * x * J1(x) * ln(x / c)`
    ///
    /// The series of the functions of the second kind and of the logarithm
    /// only converge for `|d| < c`, and their coefficients (about `c^-k`)
    /// nearly cancel, leaving much smaller ones (about `1 / k!`), so this
    /// needs a high precision.
    fn anchored_taylor(self, c: &rug::Float, n: usize) -> Vec<rug::Float> {
        assert!(matches!(self, Self::Y0 | Self::Y1));
        let pi = rug::Float::with_val(PREC, rug::float::Constant::Pi);
        let frac_2_pi = 2u32 / pi;

        let y = self.taylor(c, n);
        let j = self.first_kind().taylor(c, n);
        // ln(1 + d / c) = sum((-1)^(i + 1) * (d / c)^i / i)
        let mut ln = vec![rug::Float::new(PREC)];
        let mut ci = rug::Float::with_val(PREC, 1);
        for i in 1..n {
            ci *= c;
            let term = rug::Float::with_val(PREC, i as u32 * &ci).recip();
            ln.push(if i % 2 == 0 { -term } else { term });
        }

        // h = y - 2/π * j * ln
        let h: Vec<rug::Float> = (0..n)
            .map(|k| {
                let mut s = rug::Float::new(PREC);
                for i in 0..k {
                    s += rug::Float::with_val(PREC, &j[i] * &ln[k - i]);
                }
                y[k].clone() - s * &frac_2_pi
            })
            .collect();

        match self {
            Self::Y0 => h,
            // (c + d) * h
            _ => (0..n)
                .map(|k| {
                    let mut e = rug::Float::with_val(PREC, &h[k] * c);
                    if k >= 1 {
                        e += &h[k - 1];
                    }
                    e
                })
                .collect(),
        }
    }
}

/// Returns `value` as a Julia expression.
fn julia_float(value: &rug::Float) -> String {
    format!(
        "BigFloat(\"{}\")",
        value.to_string_radix(10, Some(JULIA_DIGITS))
    )
}

/// Returns `values` as a Julia array expression.
fn julia_array(values: &[rug::Float]) -> String {
    let items: Vec<String> = values.iter().map(julia_float).collect();
    format!("BigFloat[{}]", items.join(", "))
}

fn render_value(fkind: FloatKind, value: &rug::Float, out: &mut String) {
    render_const_value(fkind, value, out);
    out.push_str(", // ");
    render_const_dec_value(fkind, value, out);
    out.push('\n');
}

/// Returns `2^e` (as `f64`, exactly when it is normal).
fn exp2i(e: i32) -> f64 {
    (rug::Float::with_val(53, 1) << e).to_f64()
}

/// Returns the value closest to `z` such that `x - c` is exact for every `x`
/// in `[a, b]`: `c` is a multiple of `2^g`, with the smallest `g` that makes
/// all those differences fit in 53 bits.
fn exact_center(z: &rug::Float, a: f64, b: f64) -> Result<f64, String> {
    let exp_of = |v: f64| rug::Float::with_val(53, v).get_exp().unwrap() - 1;
    let g_min = exp_of(z.to_f64()) - 52;
    for g in g_min..g_min + 8 {
        let c = (rug::Float::with_val(PREC, z >> g).round() << g).to_f64();
        // In each binade `[2^e, 2^(e + 1)]`, `x` is a multiple of `2^(e - 52)`
        let exact = (exp_of(a)..=exp_of(b)).all(|e| {
            let lo = a.max(exp2i(e));
            let hi = b.min(exp2i(e + 1));
            let unit_exp = (e - 52).min(g);
            (lo - c).abs().max((hi - c).abs()) < exp2i(unit_exp + 53)
        });
        if exact {
            return Ok(c);
        }
    }
    Err(format!("no exact center found for the zero {z:e}"))
}

/// Generates a table of approximations of a Bessel function (`j0`, `j1`, `y0`
/// or `y1`) around its first `num_cells` positive zeros, as a static array of
/// `ZeroPoly` named `{FUNC}_TBL`.
///
/// The `k`-th zero (counting from 1) is close to `z = (k + offset) * π` (see
/// `Func::zero_offset`), and it is the only one in the interval
/// `[z - π/2, z + π/2)`, which is the `k`-th interval of the table, except
/// that the first one starts at `x_min`. For each interval,
///
/// `f(c + d) ~= K0 + K1 * d + ... + Kn * d^n`
///
/// where `c` is a value very close to the zero such that `d = x - c` is exact
/// in the interval, `K0 = f(c)` and the other coefficients are fitted to
/// minimize the relative error, with the first `num_hi` coefficients as
/// `F64x2` (field `k_hi`) and the next `num_lo` ones as `f64` (field `k_lo`).
///
/// The functions of the second kind have a singularity at zero, which is too
/// close to their first zeros, so `f` is the function without the
/// singularity (see `Func::anchored_taylor`) in the first `num_anchored`
/// intervals. Their fits minimize the error relative to the Bessel function.
///
/// Arguments: `func num_cells num_anchored x_min num_hi num_lo`
pub(super) fn gen_bessel_table(args: &[&str]) -> Result<String, String> {
    let [func, num_cells, num_anchored, x_min, num_hi, num_lo] = args else {
        return Err(format!("expected 6 arguments, found {}", args.len()));
    };
    let func: Func = func.parse()?;
    let (num_cells, num_anchored): (u32, u32) =
        arg_utils::parse_2_args(&[num_cells, num_anchored])?;
    let (x_min, num_hi, num_lo): (f64, usize, usize) =
        arg_utils::parse_3_args(&[x_min, num_hi, num_lo])?;
    if num_cells == 0 || num_anchored > num_cells {
        return Err(format!(
            "invalid number of intervals: {num_cells} ({num_anchored} anchored)"
        ));
    }
    if num_anchored != 0 && func == func.first_kind() {
        return Err(format!("{} has no singularity", func.name()));
    }
    if num_hi == 0 || num_lo == 0 {
        return Err("expected at least one coefficient of each kind".into());
    }

    let name = func.name();
    let poly_deg = num_hi + num_lo - 1;
    let pi = std::f64::consts::PI;
    let off = func.zero_offset();

    let mut out = String::new();
    writeln!(
        out,
        "// {}_TBL[i] has the coefficients `K0, ..., K{poly_deg}` such that",
        name.to_uppercase()
    )
    .unwrap();
    writeln!(
        out,
        "// f(c + d) ~= K0 + K1 * d + ... + K{poly_deg} * d^{poly_deg}, where `c` is very close to the"
    )
    .unwrap();
    writeln!(
        out,
        "// (i + 1)-th positive zero of {name}, for the values of the interval of that"
    )
    .unwrap();
    // "i + a" or "i - a"
    let i_plus = |a: f64| {
        if a < 0.0 {
            format!("i - {}", -a)
        } else {
            format!("i + {a}")
        }
    };
    writeln!(
        out,
        "// zero, [({}) * π, ({}) * π) (the first one starts at {x_min})",
        i_plus(off + 0.5),
        i_plus(off + 1.5),
    )
    .unwrap();
    if num_anchored == 0 {
        writeln!(out, "// f(x) = {name}(x)").unwrap();
    } else {
        let first = func.first_kind().name();
        let anchored = match func {
            Func::Y0 => format!("{name}(x) - 2/π * {first}(x) * ln(x / c)"),
            _ => format!("x * {name}(x) - 2/π * x * {first}(x) * ln(x / c)"),
        };
        writeln!(
            out,
            "// f(x) = {anchored} in the first {num_anchored} intervals"
        )
        .unwrap();
        writeln!(out, "// f(x) = {name}(x) in the other ones").unwrap();
    }
    writeln!(
        out,
        "static {}_TBL: [ZeroPoly; {num_cells}] = [",
        name.to_uppercase()
    )
    .unwrap();

    for k in 1..=num_cells {
        let a = if k == 1 {
            x_min
        } else {
            (f64::from(k) + off - 0.5) * pi
        };
        let b = (f64::from(k) + off + 0.5) * pi;
        eprintln!("{name}: interval [{a}, {b})");

        let z = func.zero(k)?;
        if !(z > a && z < b) {
            return Err(format!("the zero {z:e} is not in [{a}, {b})"));
        }
        // The interval of a value is selected with a rounded calculation, so
        // the values next to its limits can use any of the two polynomials
        // (except at `x_min`, which is checked exactly)
        let margin = 1e-9;
        let a_ext = if k == 1 { a } else { a - margin };
        let c = exact_center(&z, a_ext, b + margin)?;
        let cf = rug::Float::with_val(PREC, c);

        let anchored = k <= num_anchored;
        let taylor = if anchored {
            func.anchored_taylor(&cf, TAYLOR_TERMS)
        } else {
            func.taylor(&cf, TAYLOR_TERMS)
        };
        // The last terms must be negligible in the whole interval
        let d_max = (a - c).abs().max((b - c).abs());
        let last = rug::Float::with_val(PREC, taylor.last().unwrap() / &taylor[1])
            .abs()
            .to_f64();
        let tail = rug::Float::with_val(PREC, d_max).pow(TAYLOR_TERMS as u32 - 2) * last;
        if tail > exp2i(-140) {
            return Err(format!("not enough Taylor terms for [{a}, {b})"));
        }

        // Fit q(d) = (f(c + d) - f(c)) / d, so `K0 = f(c)` is kept exact, which
        // is what makes the relative error small close to the zero.
        //
        // The error is weighted relative to `F(c + d) - F(c)`, where `F` is the
        // Bessel function. `F(c)` is tiny, but the error relative to
        // `F(c + d)` would be unbounded at the actual zero, which is not
        // exactly `c`.
        let mut setup = format!("C = {};\n", julia_array(&taylor[1..]));
        let wfunc = if anchored {
            let first = func.first_kind().taylor(&cf, TAYLOR_TERMS);
            writeln!(setup, "J = {};", julia_array(&first)).unwrap();
            writeln!(setup, "c0 = {};", julia_float(&cf)).unwrap();
            writeln!(setup, "f0 = {};", julia_float(&taylor[0])).unwrap();
            setup.push_str("k = 2 / BigFloat(pi);\n");
            setup.push_str("lq = (x) -> iszero(x) ? 1 / c0 : log1p(x / c0) / x;\n");
            match func {
                // (y0(c + d) - y0(c)) / d = q(d) + 2/π * j0(c + d) * ln(1 + d / c) / d
                Func::Y0 => "1 / abs(fx + k * evalpoly(x, J) * lq(x))",
                // y1(c + d) = f(c + d) / (c + d) + 2/π * j1(c + d) * ln(1 + d / c),
                // so the error of `f` is divided by `c + d`, and y1(c) = f(c) / c:
                //
                // (c + d) * (y1(c + d) - y1(c)) / d
                //     = q(d) - f(c) / c + (c + d) * 2/π * j1(c + d) * ln(1 + d / c) / d
                _ => "1 / abs((c0 * fx - f0) / c0 + (c0 + x) * k * evalpoly(x, J) * lq(x))",
            }
        } else {
            // `f` is the Bessel function, so (f(c + d) - f(c)) / d = q(d)
            "1 / abs(fx)"
        };
        let coeffs = julia::run_remez_with_setup(
            &setup,
            "evalpoly(x, C)",
            wfunc,
            (a - c, b - c),
            poly_deg as i32 - 1,
        );

        writeln!(out, "// [{a}, {b})").unwrap();
        writeln!(out, "ZeroPoly {{").unwrap();
        out.push_str("c: ");
        render_value(FloatKind::F64, &cf, &mut out);
        let all = std::iter::once(&taylor[0]).chain(coeffs.iter());
        writeln!(out, "k_hi: [").unwrap();
        for coeff in all.clone().take(num_hi) {
            render_value(FloatKind::F64x2, coeff, &mut out);
        }
        writeln!(out, "],").unwrap();
        writeln!(out, "k_lo: [").unwrap();
        for coeff in all.skip(num_hi) {
            render_value(FloatKind::F64, coeff, &mut out);
        }
        writeln!(out, "],").unwrap();
        writeln!(out, "}},").unwrap();
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

/// Generates an approximation of a Bessel function for `0 < x <= x_max`, with
/// a polynomial `P(u)` of `u = x^2`:
///
/// * `j0(x) ~= 1 + u * P(u)`, with `P(u) = K1 + K2 * u + ...`
/// * `j1(x) ~= x * P(u)`, with `P(u) = K0 + K1 * u + ...`
/// * `y0(x) ~= 2/π * ln(x) * j0(x) + P(u)`, with `P(u) = K0 + K1 * u + ...`
/// * `y1(x) ~= 2/π * ln(x) * j1(x) - 2 / (π * x) + x * P(u)`, with
///   `P(u) = K0 + K1 * u + ...`
///
/// minimizing the error relative to the Bessel function.
///
/// Arguments: `func fkinds num_coeffs x_max`
pub(super) fn gen_bessel_small_poly(args: &[&str]) -> Result<String, String> {
    let (func, fkinds, num_coeffs, x_max): (Func, CoeffKinds, i32, f64) =
        arg_utils::parse_4_args(args)?;
    if num_coeffs < 1 {
        return Err("expected at least 1 coefficient".into());
    }

    let mut out = String::new();

    // The fitted functions of `u = x^2` are not defined at 0, so the range
    // starts at a tiny value
    let range = (1e-40, x_max * x_max);
    let (fit, wfunc, first) = match func {
        Func::J0 => (
            "(SpecialFunctions.besselj0(sqrt(x)) - 1) / x",
            "x / SpecialFunctions.besselj0(sqrt(x))",
            1,
        ),
        Func::J1 => (
            "begin s = sqrt(x); SpecialFunctions.besselj1(s) / s end",
            "1 / fx",
            0,
        ),
        Func::Y0 => (
            "begin \
            s = sqrt(x); \
            SpecialFunctions.bessely0(s) - 2 / BigFloat(pi) * log(s) * SpecialFunctions.besselj0(s) \
            end",
            "1 / abs(SpecialFunctions.bessely0(sqrt(x)))",
            0,
        ),
        Func::Y1 => (
            "begin \
            s = sqrt(x); \
            (SpecialFunctions.bessely1(s) - 2 / BigFloat(pi) * log(s) * SpecialFunctions.besselj1(s) \
            + 2 / (BigFloat(pi) * s)) / s \
            end",
            "begin s = sqrt(x); s / abs(SpecialFunctions.bessely1(s)) end",
            0,
        ),
    };

    julia::run_and_render_remez(
        &fkinds,
        fit,
        wfunc,
        range,
        num_coeffs - 1,
        first,
        "K",
        &mut out,
    );

    Ok(out)
}

/// Asymptotic series of the modulus and the phase of the Bessel functions of
/// order `nu`, for large `x`, where `s = 1 / x^2`:
///
/// * `J(x)^2 + Y(x)^2 ~ 2 / (π * x) * sum(modulus_sq[k] * s^k)`
/// * `J(x) = M(x) * cos(θ(x))` and `Y(x) = M(x) * sin(θ(x))`, with
///   `θ(x) ~ x - (nu / 2 + 1/4) * π + sum(phase[k] * s^k) / x`
struct AsymptSeries {
    modulus_sq: Vec<rug::Float>,
    phase: Vec<rug::Float>,
}

impl AsymptSeries {
    /// Number of terms of the series.
    const TERMS: usize = 45;

    /// Returns the series for the order `nu`, to be used with `x >= x0`.
    ///
    /// The terms of the series decrease up to about the `x`-th one and then
    /// they grow, so the series are only accurate with large values: this
    /// fails if the last term of the phase with `x0` is not negligible
    /// (2^-140, which is much less than the errors of the polynomials fitted
    /// to these series).
    fn new(nu: u32, x0: f64) -> Result<Self, String> {
        if nu > 1 {
            return Err(format!("invalid order: {nu}"));
        }
        if !(x0 >= 1.0 && x0.is_finite()) {
            return Err(format!("invalid x0: {x0}"));
        }
        let series = Self::with_order(nu);
        let last: rug::Float = series.phase[Self::TERMS - 1].clone().abs()
            / rug::Float::with_val(PREC, x0).pow(2 * Self::TERMS as u32 - 1);
        if last > rug::Float::with_val(PREC, 1) >> 140u32 {
            return Err(format!("the asymptotic series do not converge with {x0}"));
        }
        Ok(series)
    }

    fn with_order(nu: u32) -> Self {
        // modulus_sq[k] = modulus_sq[k-1] * (2k - 1) / (2k) * (μ - (2k - 1)^2) / 4,
        // with μ = 4 * ν^2 (the coefficients are rational, so the precision
        // used here makes them exact for all practical purposes)
        let mu = i64::from(4 * nu * nu);
        let mut modulus_sq = vec![rug::Float::with_val(PREC, 1)];
        for k in 1..=Self::TERMS as i64 {
            let prev = modulus_sq.last().unwrap();
            let next = rug::Float::with_val(PREC, prev * (2 * k - 1)) / (2 * k)
                * (mu - (2 * k - 1).pow(2))
                / 4i64;
            modulus_sq.push(next);
        }

        // θ'(x) = 2 / (π * x * M(x)^2) = 1 / sum(modulus_sq[k] * s^k) = sum(inv[k] * s^k)
        let mut inv = vec![rug::Float::with_val(PREC, 1)];
        for k in 1..=Self::TERMS {
            let mut sum = rug::Float::new(PREC);
            for j in 1..=k {
                sum += rug::Float::with_val(PREC, &modulus_sq[j] * &inv[k - j]);
            }
            inv.push(-sum);
        }

        // Integrating the terms after the first one:
        // phase[k] = -inv[k + 1] / (2 * k + 1)
        let phase = (0..Self::TERMS)
            .map(|k| -inv[k + 1].clone() / (2 * k as u32 + 1))
            .collect();

        modulus_sq.truncate(Self::TERMS);
        Self { modulus_sq, phase }
    }
}

/// Fits `sum(phase[k] * s^k)` (see `AsymptSeries`) for
/// `0 <= s <= 1 / x0^2` with a polynomial with `num_coeffs` coefficients,
/// where the first one is the one of the series and the others minimize the
/// absolute error of the phase.
fn fit_phase(nu: u32, num_coeffs: i32, x0: f64) -> Result<Vec<rug::Float>, String> {
    let series = AsymptSeries::new(nu, x0)?;
    // sum(phase[k] * s^k) = phase[0] + s * R(s), and the error of the phase is
    // the one of `R(s)` times `s / x = s^(3/2)`
    let setup = format!("C = {};", julia_array(&series.phase[1..]));
    let range = (1e-40, 1.0 / (x0 * x0));
    let fitted = julia::run_remez_with_setup(
        &setup,
        "evalpoly(x, C)",
        "x * sqrt(x)",
        range,
        num_coeffs - 2,
    );
    Ok(std::iter::once(series.phase[0].clone())
        .chain(fitted)
        .collect())
}

/// Generates an approximation of the phase of the Bessel functions of order
/// `nu` (0 or 1), for `x >= x0`:
///
/// `θ(x) ~= x - (nu / 2 + 1/4) * π + (K0 + K1 / x^2 + ... + Kn / x^(2 * n)) / x`
///
/// where `K0` is the coefficient of the asymptotic series and `K1, K2, ...`
/// are fitted to minimize the absolute error.
///
/// Arguments: `nu fkinds num_coeffs x0`, where `num_coeffs` includes `K0`.
pub(super) fn gen_bessel_phase_poly(args: &[&str]) -> Result<String, String> {
    let (nu, fkinds, num_coeffs, x0): (u32, CoeffKinds, i32, f64) = arg_utils::parse_4_args(args)?;
    if num_coeffs < 2 {
        return Err("expected at least 2 coefficients".into());
    }

    let mut out = String::new();
    for (i, coeff) in fit_phase(nu, num_coeffs, x0)?.into_iter().enumerate() {
        render_const(fkinds.get(i), &format!("K{i}"), coeff, &mut out);
    }

    Ok(out)
}

/// Generates the coefficients of `gen_bessel_phase_poly` as a static array of
/// fixed point numbers named `PHASE{nu}_FIXED`, to calculate the phase with
/// integer arithmetic, in units of π/128:
///
/// `θ(x) * 128/π ~= x * 128/π - (nu / 2 + 1/4) * 128 + sum(D[k] * σ^k) / x`
///
/// with `σ = 2^10 / x^2` and `D[k] = 128/π * K[k] * 2^(-10 * k)` (which
/// decrease quickly), in units of `2^-frac_bits`.
///
/// Arguments: `nu num_coeffs x0 frac_bits`
pub(super) fn gen_bessel_phase_fixed(args: &[&str]) -> Result<String, String> {
    let (nu, num_coeffs, x0, frac_bits): (u32, i32, f64, u32) = arg_utils::parse_4_args(args)?;
    if num_coeffs < 2 {
        return Err("expected at least 2 coefficients".into());
    }
    if x0 < 32.0 {
        // σ must be less than one
        return Err(format!("invalid x0: {x0}"));
    }
    // It also checks the order
    let coeffs = fit_phase(nu, num_coeffs, x0)?;

    let mut out = String::new();

    writeln!(
        out,
        "// PHASE{nu}_FIXED[k] = 128/π * K[k] * 2^(-10 * k), in units of 2^-{frac_bits}, where"
    )
    .unwrap();
    writeln!(
        out,
        "// θ(x) ~= x - {} + (K[0] + K[1] / x^2 + ... + K[{}] / x^{}) / x, for x >= {x0}",
        ["π/4", "3 * π/4"][nu as usize],
        num_coeffs - 1,
        2 * (num_coeffs - 1),
    )
    .unwrap();
    writeln!(out, "static PHASE{nu}_FIXED: [i128; {num_coeffs}] = [").unwrap();

    let scale = 128u32 / rug::Float::with_val(PREC, rug::float::Constant::Pi);
    for (k, coeff) in coeffs.into_iter().enumerate() {
        let d = (coeff * &scale) >> (10 * k as u32);
        let fixed = rug::Float::with_val(PREC, &d << frac_bits)
            .round()
            .to_integer()
            .unwrap();
        let Some(fixed) = fixed.to_i128() else {
            return Err(format!("coefficient {k} does not fit in 128 bits"));
        };
        let sign = if fixed < 0 { "-" } else { "" };
        writeln!(
            out,
            "    {sign}0x{:032X}, // {:e}",
            fixed.unsigned_abs(),
            d.to_f64()
        )
        .unwrap();
    }
    out.push_str("];\n");

    Ok(out)
}

/// Generates an approximation of the modulus of the Bessel functions of order
/// `nu` (0 or 1), for `x >= x0`:
///
/// `M(x) ~= sqrt(2 / (π * x)) * (1 + K1 / x^2 + ... + Kn / x^(2 * n))`
///
/// where `K1, K2, ...` are fitted to minimize the relative error.
///
/// Arguments: `nu fkinds num_coeffs x0`
pub(super) fn gen_bessel_modulus_poly(args: &[&str]) -> Result<String, String> {
    let (nu, fkinds, num_coeffs, x0): (u32, CoeffKinds, i32, f64) = arg_utils::parse_4_args(args)?;
    if num_coeffs < 1 {
        return Err("expected at least 1 coefficient".into());
    }

    let mut out = String::new();

    // sqrt(sum(modulus_sq[k] * s^k)) = 1 + s * R(s), and the error is the one
    // of `R(s)` times `s`
    let series = AsymptSeries::new(nu, x0)?;
    let setup = format!("C = {};", julia_array(&series.modulus_sq));
    let range = (1e-40, 1.0 / (x0 * x0));
    let fitted = julia::run_remez_with_setup(
        &setup,
        "(sqrt(evalpoly(x, C)) - 1) / x",
        "x",
        range,
        num_coeffs - 1,
    );
    for (i, coeff) in fitted.into_iter().enumerate() {
        render_const(fkinds.get(i), &format!("K{}", i + 1), coeff, &mut out);
    }

    Ok(out)
}

/// Generates the coefficients of the Debye expansions of the Bessel functions
/// of large order `nu`, up to the term of degree `num_terms`.
///
/// The expansions have the factor `sum(U[k](p) / nu^k)` (with `U[0] = 1`),
/// whose logarithm is `sum(W[k](p) / nu^k)` (with `k >= 1`), where `W[k]` has
/// the powers `p^k, p^(k + 2), ..., p^(3 * k)`. With `a = p^3 / nu` and
/// `b = p / nu`, it is
///
/// `P(a, b) = sum(C[k][i] * a^(k - i) * b^i)`, for `k >= 1` and `0 <= i <= k`
///
/// where `C[k][i]` is the coefficient of `p^(3 * k - 2 * i)` in `W[k]`. The
/// terms of degree 1 are `DEBYE_A1 * a + DEBYE_B1 * b`, and `DEBYE_POLY` has
/// the coefficients of the other ones (`C[m + i][i]` for `a^m * b^i`),
/// grouped by `i`, with zeros for the terms of degree less than 2.
///
/// Arguments: `num_terms`
pub(super) fn gen_bessel_debye_poly(args: &[&str]) -> Result<String, String> {
    let [num_terms] = args else {
        return Err(format!("expected 1 argument, found {}", args.len()));
    };
    let num_terms: u32 = num_terms
        .parse()
        .map_err(|e| format!("failed to parse number of terms {num_terms:?}: {e}"))?;
    if num_terms < 2 {
        return Err(format!("expected at least 2 terms, found {num_terms}"));
    }
    let k_max = num_terms as usize;
    let zero = || rug::Float::new(PREC);
    let u = debye_u(k_max);

    // The logarithm, from S' / S, where S = sum(U[k] * e^k):
    // k * W[k] = k * U[k] - sum(j * W[j] * U[k - j]), for 1 <= j < k
    let mut w: Vec<Vec<rug::Float>> = vec![Vec::new()];
    for k in 1..=k_max {
        let mut acc: Vec<rug::Float> = u[k]
            .iter()
            .map(|c| rug::Float::with_val(PREC, c * k as u32))
            .collect();
        for j in 1..k {
            for (ia, ca) in w[j].iter().enumerate() {
                for (ib, cb) in u[k - j].iter().enumerate() {
                    acc[ia + ib] -= rug::Float::with_val(PREC, ca * cb) * j as u32;
                }
            }
        }
        w.push(acc.into_iter().map(|c| c / k as u32).collect());
    }

    // C[k][i], checking the powers of `W[k]`
    let tiny = rug::Float::with_val(PREC, 1) >> (PREC / 2);
    let coeff = |k: usize, i: usize| -> Result<&rug::Float, String> {
        for (j, c) in w[k].iter().enumerate() {
            let expected = j >= k && (j - k) % 2 == 0;
            if !expected && c.clone().abs() > tiny {
                return Err(format!("unexpected power p^{j} in W[{k}]"));
            }
        }
        Ok(&w[k][3 * k - 2 * i])
    };

    let mut out = String::new();
    writeln!(
        out,
        "// P(a, b) = sum(C[k][i] * a^(k - i) * b^i), for 1 <= k <= {k_max} and 0 <= i <= k"
    )
    .unwrap();
    render_const(FloatKind::F64x2, "DEBYE_A1", coeff(1, 0)?.clone(), &mut out);
    render_const(FloatKind::F64x2, "DEBYE_B1", coeff(1, 1)?.clone(), &mut out);

    let len = (k_max + 1) * (k_max + 2) / 2;
    writeln!(
        out,
        "// DEBYE_POLY has C[m + i][i] for a^m * b^i, with m + i >= 2, grouped by i"
    )
    .unwrap();
    writeln!(out, "static DEBYE_POLY: [f64; {len}] = [").unwrap();
    for i in 0..=k_max {
        writeln!(out, "    // b^{i}").unwrap();
        for m in 0..=(k_max - i) {
            let k = m + i;
            let c = if k < 2 { zero() } else { coeff(k, i)?.clone() };
            out.push_str("    ");
            render_value(FloatKind::F64, &c, &mut out);
        }
    }
    out.push_str("];\n");

    Ok(out)
}

/// Returns the coefficients of the powers of `p` of the polynomials `U[k](p)`
/// of the Debye expansions (DLMF 10.41.10), for `0 <= k <= k_max`.
fn debye_u(k_max: usize) -> Vec<Vec<rug::Float>> {
    // U[k + 1](p) = p^2 * (1 - p^2) / 2 * U[k]'(p)
    //     + 1/8 * integral(0..p, (1 - 5 * t^2) * U[k](t) dt)
    let mut u = vec![vec![rug::Float::with_val(PREC, 1)]];
    for k in 0..k_max {
        let prev = &u[k];
        let mut next = vec![rug::Float::new(PREC); prev.len() + 3];
        for (j, c) in prev.iter().enumerate() {
            let jf = j as u32;
            // j * c * p^(j - 1) * (p^2 - p^4) / 2
            let t = rug::Float::with_val(PREC, c * jf) / 2u32;
            next[j + 1] += &t;
            next[j + 3] -= &t;
            // c * (p^(j + 1) / (j + 1) - 5 * p^(j + 3) / (j + 3)) / 8
            next[j + 1] += rug::Float::with_val(PREC, c / (8 * (jf + 1)));
            next[j + 3] -= rug::Float::with_val(PREC, c * 5u32) / (8 * (jf + 3));
        }
        u.push(next);
    }
    u
}

/// Generates a table of the Airy functions and their derivatives,
/// `[Ai(t), Ai'(t), Bi(t), Bi'(t)]` for `t = t_min + i / inv_step`, up to
/// `t_max`.
///
/// Arguments: `t_min t_max inv_step`
pub(super) fn gen_bessel_airy_table(args: &[&str]) -> Result<String, String> {
    let (t_min, t_max, inv_step): (i32, i32, u32) = arg_utils::parse_3_args(args)?;
    if t_min >= t_max || inv_step == 0 {
        return Err(format!("invalid range: {t_min} {t_max} {inv_step}"));
    }
    let num = (t_max - t_min) as u32 * inv_step + 1;

    // Ai(0) = 3^(-2/3) / Γ(2/3), Ai'(0) = -3^(-1/3) / Γ(1/3),
    // Bi(0) = sqrt(3) * Ai(0) and Bi'(0) = -sqrt(3) * Ai'(0)
    let third = rug::Float::with_val(PREC, 1) / 3u32;
    let three = rug::Float::with_val(PREC, 3);
    let two_thirds = rug::Float::with_val(PREC, &third * 2u32);
    let ai0 = rug::Float::with_val(PREC, (&three).pow(-two_thirds.clone())) / two_thirds.gamma();
    let dai0 = -rug::Float::with_val(PREC, (&three).pow(-third.clone()))
        / rug::Float::with_val(PREC, third.gamma_ref());
    let sqrt3 = three.sqrt();
    let bi0 = rug::Float::with_val(PREC, &ai0 * &sqrt3);
    let dbi0 = -rug::Float::with_val(PREC, &dai0 * &sqrt3);

    // The Wronskian below checks the scale of `Bi`, but not a multiple of
    // `Ai` added to it (which does not change it). That would change the
    // phase of its oscillations with `t < 0`, so both functions are also
    // checked with their asymptotic expansions at `t = -26`, whose error is
    // about 2^-255 there.
    let t_check = rug::Float::with_val(PREC, -26);
    let (ai_asym, bi_asym) = airy_neg_asymptotic(&rug::Float::with_val(PREC, -&t_check));
    let (ai_check, _) = airy_maclaurin(&ai0, &dai0, &t_check);
    let (bi_check, _) = airy_maclaurin(&bi0, &dbi0, &t_check);
    let modulus = rug::Float::with_val(PREC, ai_check.hypot_ref(&bi_check));
    let tol = rug::Float::with_val(PREC, &modulus >> 240u32);
    if (ai_check - ai_asym).abs() > tol || (bi_check - bi_asym).abs() > tol {
        return Err(format!("inaccurate Airy functions at {}", t_check.to_f64()));
    }

    let mut out = String::new();
    writeln!(
        out,
        "// AIRY_TBL[i] = [Ai(t), Ai'(t), Bi(t), Bi'(t)], with t = {t_min} + i / {inv_step}"
    )
    .unwrap();
    writeln!(out, "static AIRY_TBL: [[F64x2; 4]; {num}] = [").unwrap();
    for i in 0..num {
        let t = rug::Float::with_val(PREC, t_min) + rug::Float::with_val(PREC, i) / inv_step;
        let (ai, dai) = airy_maclaurin(&ai0, &dai0, &t);
        let (bi, dbi) = airy_maclaurin(&bi0, &dbi0, &t);
        // Check the series with MPFR (which only has `Ai`), and the others
        // with their Wronskian, Ai(t) * Bi'(t) - Ai'(t) * Bi(t) = 1 / π
        let tol = rug::Float::with_val(PREC, 1) >> 1000u32;
        let ai_ref = t.clone().ai();
        let diff = rug::Float::with_val(PREC, &ai - &ai_ref) / &ai_ref;
        let w = rug::Float::with_val(PREC, &ai * &dbi) - rug::Float::with_val(PREC, &dai * &bi);
        let w_diff = w * rug::Float::with_val(PREC, rug::float::Constant::Pi) - 1u32;
        if diff.abs() > tol || w_diff.abs() > tol {
            return Err(format!("inaccurate Airy functions at {}", t.to_f64()));
        }
        writeln!(out, "    // t = {}", t.to_f64()).unwrap();
        out.push_str("    [\n");
        for v in [&ai, &dai, &bi, &dbi] {
            out.push_str("        ");
            render_value(FloatKind::F64x2, v, &mut out);
        }
        out.push_str("    ],\n");
    }
    out.push_str("];\n");

    Ok(out)
}

/// Returns the coefficients `u[k]` of the asymptotic expansions of the Airy
/// functions (DLMF 9.7.2), for `0 <= k <= k_max`.
fn airy_u(k_max: usize) -> Vec<rug::Float> {
    // u[k] = (6k - 5) * (6k - 3) * (6k - 1) / ((2k - 1) * 216 * k) * u[k - 1]
    let mut u = vec![rug::Float::with_val(PREC, 1)];
    for k in 1..=k_max as u32 {
        let next = rug::Float::with_val(PREC, &u[k as usize - 1])
            * (6 * k - 5)
            * (6 * k - 3)
            * (6 * k - 1)
            / ((2 * k - 1) * 216 * k);
        u.push(next);
    }
    u
}

/// Returns `(Ai(-x), Bi(-x))` from their asymptotic expansions (DLMF 9.7.9
/// and 9.7.11), truncated before their smallest term, for large `x > 0`.
///
/// With `ζ = 2/3 * x^(3/2)`, they are
///
/// * `Ai(-x) ~= (cos(ζ - π/4) * P + sin(ζ - π/4) * Q) / (sqrt(π) * x^(1/4))`
/// * `Bi(-x) ~= (cos(ζ - π/4) * Q - sin(ζ - π/4) * P) / (sqrt(π) * x^(1/4))`
///
/// where `P = sum((-1)^k * u[2k] / ζ^(2k))` and
/// `Q = sum((-1)^k * u[2k + 1] / ζ^(2k + 1))`.
fn airy_neg_asymptotic(x: &rug::Float) -> (rug::Float, rug::Float) {
    let zeta = rug::Float::with_val(PREC, x.sqrt_ref()) * x * 2u32 / 3u32;
    // The terms decrease up to `k` about `2 * ζ`
    let k_max = (zeta.to_f64() * 2.0) as usize;
    let u = airy_u(k_max + 1);
    let (mut p, mut q) = (rug::Float::new(PREC), rug::Float::new(PREC));
    let mut zeta_pow = rug::Float::with_val(PREC, 1);
    let mut prev = None;
    for (k, uk) in u.iter().enumerate() {
        let term = rug::Float::with_val(PREC, uk / &zeta_pow);
        if prev.is_some_and(|prev: rug::Float| term.clone().abs() >= prev) {
            break;
        }
        match k % 4 {
            0 => p += &term,
            1 => q += &term,
            2 => p -= &term,
            _ => q -= &term,
        }
        prev = Some(term.abs());
        zeta_pow *= &zeta;
    }

    let pi = rug::Float::with_val(PREC, rug::float::Constant::Pi);
    let phase = rug::Float::with_val(PREC, &zeta - rug::Float::with_val(PREC, &pi / 4u32));
    let (sin, cos) = phase.sin_cos(rug::Float::new(PREC));
    let scale = rug::Float::with_val(PREC, pi.sqrt() * x.clone().sqrt().sqrt()).recip();
    let ai =
        (rug::Float::with_val(PREC, &cos * &p) + rug::Float::with_val(PREC, &sin * &q)) * &scale;
    let bi =
        (rug::Float::with_val(PREC, &cos * &q) - rug::Float::with_val(PREC, &sin * &p)) * &scale;
    (ai, bi)
}

/// Returns `(y(t), y'(t))` from the Maclaurin series of the solution of
/// `y'' = t * y` with `y(0) = y0` and `y'(0) = y1`.
///
/// Its terms are up to about `exp(2/3 * |t|^(3/2))` times larger than the
/// Airy functions, which are as small as the inverse of that (Ai with
/// `t > 0`), so the precision must be much higher than the one of the
/// results.
fn airy_maclaurin(y0: &rug::Float, y1: &rug::Float, t: &rug::Float) -> (rug::Float, rug::Float) {
    // y = sum(c[m] * t^m), with (m + 2) * (m + 1) * c[m + 2] = c[m - 1]
    let tiny = rug::Float::with_val(PREC, 1) >> (PREC + 64);
    let mut c = vec![y0.clone(), y1.clone(), rug::Float::new(PREC)];
    let mut y = rug::Float::with_val(PREC, y1 * t) + y0;
    let mut dy = y1.clone();
    // t^(m - 1)
    let mut t_pow = t.clone();
    let mut small = 0;
    let mut m = 2;
    while small < 3 {
        if m >= 3 {
            let cm = rug::Float::with_val(PREC, &c[m - 3] / ((m as u32 - 1) * m as u32));
            c.push(cm);
        }
        let d = rug::Float::with_val(PREC, &c[m] * &t_pow) * m as u32;
        t_pow *= t;
        let term = rug::Float::with_val(PREC, &c[m] * &t_pow);
        if term.clone().abs() < tiny && d.clone().abs() < tiny {
            small += 1;
        } else {
            small = 0;
        }
        y += term;
        dy += d;
        m += 1;
    }
    (y, dy)
}

/// Generates an approximation of the coefficient `A[k](ζ)` (`kind = A`) or
/// `B[k](ζ)` (`kind = B`) of Olver's uniform asymptotic expansions of the
/// Bessel functions of large order (DLMF 10.20.4), for
/// `zeta_min <= ζ <= zeta_max`, as a polynomial of degree `degree` (from
/// Chebyshev interpolation), as the array `OLVER_{kind}{k}`.
///
/// Their definitions (DLMF 10.20.10 and 10.20.11) are, with
/// `τ = 1 - z^2`, `ρ = ζ / τ` (which is positive) and `c[m][l]` the
/// coefficient of `p^l` in `U[m](p)` (see `debye_u`):
///
/// * `A[k] = sum((3/2)^j * v[j] * c[2k - j][l] * ρ^(-3j/2) * τ^(-(3j + l)/2))`
/// * `B[k] = -sum((3/2)^j * u[j] * c[2k - j + 1][l] * ρ^(-(1 + 3j)/2) * τ^(-(1 + 3j + l)/2))`
///
/// where the powers of `τ` are integers, so they also hold for `ζ < 0`
/// (`z > 1`). Their terms nearly cancel when `ζ` is small, which the high
/// precision of the calculations makes up for.
///
/// Arguments: `kind k degree zeta_min zeta_max`
pub(super) fn gen_bessel_olver_poly(args: &[&str]) -> Result<String, String> {
    let (kind, k, degree, zeta_min, zeta_max): (String, usize, usize, f64, f64) =
        arg_utils::parse_5_args(args)?;
    let second = match kind.as_str() {
        "A" => false,
        "B" => true,
        _ => return Err(format!("invalid kind: {kind:?}")),
    };
    if zeta_min >= zeta_max {
        return Err(format!("invalid range: [{zeta_min}, {zeta_max}]"));
    }

    let u = debye_u(2 * k + 1);
    // u[j] and v[j] of the asymptotic expansions of the Airy functions
    let uj = airy_u(2 * k + 1);
    let vj: Vec<rug::Float> = (0..uj.len() as u32)
        .map(|j| match j {
            0 => rug::Float::with_val(PREC, 1),
            _ => rug::Float::with_val(PREC, -&uj[j as usize]) * (6 * j + 1) / (6 * j - 1),
        })
        .collect();
    let coeff = |zeta: &rug::Float| olver_coeff(&u, &uj, &vj, second, k, zeta);

    let poly = cheb_fit(&coeff, zeta_min, zeta_max, degree);

    // The error of the polynomial with the coefficients rounded to `f64`
    let rounded: Vec<rug::Float> = poly
        .iter()
        .map(|c| rug::Float::with_val(PREC, c.to_f64()))
        .collect();
    let mut max_err = rug::Float::new(PREC);
    for i in 0..=1000u32 {
        let zeta = rug::Float::with_val(PREC, zeta_min)
            + rug::Float::with_val(PREC, zeta_max - zeta_min) * i / 1000u32;
        if zeta.is_zero() {
            continue;
        }
        let approx = rounded
            .iter()
            .rev()
            .fold(rug::Float::new(PREC), |acc, c| acc * &zeta + c);
        let err = (approx - coeff(&zeta)).abs();
        if err > max_err {
            max_err = err;
        }
    }
    eprintln!(
        "error = {:e} = 2^({})",
        max_err.to_f64(),
        max_err.clone().log2().to_f64()
    );

    let name = format!("OLVER_{kind}{k}");
    let mut out = String::new();
    writeln!(out, "// {name}[i] is the coefficient of ζ^i").unwrap();
    writeln!(out, "const {name}: [f64; {}] = [", degree + 1).unwrap();
    for c in rounded.iter() {
        out.push_str("    ");
        render_value(FloatKind::F64, c, &mut out);
    }
    out.push_str("];\n");
    Ok(out)
}

/// Returns `A[k](ζ)` (`second = false`) or `B[k](ζ)`, for `ζ != 0` (see
/// `gen_bessel_olver_poly`).
fn olver_coeff(
    u: &[Vec<rug::Float>],
    uj: &[rug::Float],
    vj: &[rug::Float],
    second: bool,
    k: usize,
    zeta: &rug::Float,
) -> rug::Float {
    // (2/3) * |ζ|^(3/2) = atanh(y) - y with y = sqrt(1 - z^2) (ζ > 0), or
    // y - atan(y) with y = sqrt(z^2 - 1) (ζ < 0)
    let below = zeta.is_sign_positive();
    let a = rug::Float::with_val(PREC, zeta.abs_ref());
    let c = rug::Float::with_val(PREC, a.sqrt_ref()) * &a * 2u32 / 3u32;
    let f = |y: &rug::Float| {
        if below {
            rug::Float::with_val(PREC, y.atanh_ref()) - y
        } else {
            y - rug::Float::with_val(PREC, y.atan_ref())
        }
    };
    let df = |y: &rug::Float| {
        let y2 = rug::Float::with_val(PREC, y.square_ref());
        if below {
            rug::Float::with_val(PREC, &y2 / (1 - y2.clone()))
        } else {
            rug::Float::with_val(PREC, &y2 / (1 + y2.clone()))
        }
    };
    let mut y = rug::Float::with_val(PREC, &c * 3u32).cbrt();
    for _ in 0..200 {
        let step = (f(&y) - &c) / df(&y);
        y -= &step;
        if step.abs() < (rug::Float::with_val(PREC, 1) >> (PREC - 16)) {
            break;
        }
    }
    let y2 = rug::Float::with_val(PREC, y.square_ref());
    let tau = if below { y2 } else { -y2 };
    let rho_sqrt = rug::Float::with_val(PREC, zeta / &tau).sqrt();

    let mut sum = rug::Float::new(PREC);
    let j_max = if second { 2 * k + 1 } else { 2 * k };
    let three_halves = rug::Float::with_val(PREC, 3) / 2u32;
    for j in 0..=j_max {
        let (m, a, r_exp) = if second {
            (2 * k + 1 - j, &uj[j], 1 + 3 * j)
        } else {
            (2 * k - j, &vj[j], 3 * j)
        };
        let scale = rug::Float::with_val(PREC, (&three_halves).pow(j as u32))
            * a
            * rug::Float::with_val(PREC, (&rho_sqrt).pow(-(r_exp as i32)));
        for (l, cl) in u[m].iter().enumerate() {
            if cl.is_zero() {
                continue;
            }
            let e = r_exp + l;
            assert!(e % 2 == 0, "odd power of τ");
            let t = rug::Float::with_val(PREC, (&tau).pow(-((e / 2) as i32)));
            sum += t * cl * &scale;
        }
    }
    if second { -sum } else { sum }
}

/// Returns the coefficients (of the powers of `x`) of the polynomial of
/// degree `degree` that interpolates `f` at the Chebyshev nodes of
/// `[a, b]`.
fn cheb_fit(
    f: &dyn Fn(&rug::Float) -> rug::Float,
    a: f64,
    b: f64,
    degree: usize,
) -> Vec<rug::Float> {
    let nn = degree + 1;
    let pi = rug::Float::with_val(PREC, rug::float::Constant::Pi);
    let mid = rug::Float::with_val(PREC, a + b) / 2u32;
    let half = rug::Float::with_val(PREC, b - a) / 2u32;

    // Values at the nodes cos(π * (i + 1/2) / nn)
    let nodes: Vec<rug::Float> = (0..nn)
        .map(|i| (rug::Float::with_val(PREC, &pi * (2 * i as u32 + 1)) / (2 * nn as u32)).cos())
        .collect();
    let values: Vec<rug::Float> = nodes
        .iter()
        .map(|s| {
            f(&rug::Float::with_val(
                PREC,
                &mid + rug::Float::with_val(PREC, s * &half),
            ))
        })
        .collect();

    // Chebyshev coefficients
    let mut cheb = Vec::with_capacity(nn);
    for j in 0..nn {
        let mut acc = rug::Float::new(PREC);
        for (s, v) in nodes.iter().zip(&values) {
            let tj = rug::Float::with_val(PREC, s.acos_ref()) * j as u32;
            acc += tj.cos() * v;
        }
        acc = acc * 2u32 / nn as u32;
        if j == 0 {
            acc /= 2u32;
        }
        cheb.push(acc);
    }

    // T[j](s) with s = (x - mid) / half, as polynomials of `x`, from
    // T[0] = 1, T[1] = s and T[j + 1] = 2 * s * T[j] - T[j - 1]
    let s0 = rug::Float::with_val(PREC, -&mid) / &half;
    let s1 = rug::Float::with_val(PREC, half.recip_ref());
    let mut t: Vec<Vec<rug::Float>> = vec![vec![rug::Float::with_val(PREC, 1)]];
    for j in 0..degree {
        let factor = if j == 0 { 1u32 } else { 2u32 };
        let mut next = vec![rug::Float::new(PREC); j + 2];
        for (i, c) in t[j].iter().enumerate() {
            let c = rug::Float::with_val(PREC, c * factor);
            next[i] += rug::Float::with_val(PREC, &c * &s0);
            next[i + 1] += c * &s1;
        }
        if j >= 1 {
            for (r, c) in next.iter_mut().zip(&t[j - 1]) {
                *r -= c;
            }
        }
        t.push(next);
    }

    // sum(cheb[j] * T[j])
    let mut result = vec![rug::Float::new(PREC); nn];
    for (tj, cj) in t.iter().zip(&cheb) {
        for (r, c) in result.iter_mut().zip(tj) {
            *r += rug::Float::with_val(PREC, c * cj);
        }
    }
    result
}
