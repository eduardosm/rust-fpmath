pub(crate) fn gen_frac_2_pi_large() {
    let num_words = 66;
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(num_words * 24);
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.f64_div(2.0, None, dev_mpfr::Rnd::N);

    println!("const FRAC_2_PI_LARGE: &[u32] = &[");
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            print!("    ");
        } else {
            print!(" ");
        }
        tmp.mul_2ui(None, 24, dev_mpfr::Rnd::N);
        let word = tmp.get_ui(dev_mpfr::Rnd::Z);
        print!("0x{word:06X},");
        tmp.sub_ui(None, word, dev_mpfr::Rnd::N);
        if (i % words_per_line) == (words_per_line - 1) {
            println!();
        }
    }
    if (num_words % words_per_line) != 0 {
        println!();
    }
    println!("];");
}

pub(crate) fn gen_frac_pi_2_medium() {
    let num_words = 8;
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(num_words * 24);
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.mul_f64(None, 0.5, dev_mpfr::Rnd::N);

    println!("const FRAC_PI_2_MEDIUM: &[u32] = &[");
    tmp.mul_f64(None, 0.5, dev_mpfr::Rnd::N);
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            print!("    ");
        } else {
            print!(" ");
        }
        tmp.mul_2ui(None, 24, dev_mpfr::Rnd::N);
        let word = tmp.get_ui(dev_mpfr::Rnd::Z);
        print!("0x{word:06X},");
        tmp.sub_ui(None, word, dev_mpfr::Rnd::N);
        if (i % words_per_line) == (words_per_line - 1) {
            println!();
        }
    }
    if (num_words % words_per_line) != 0 {
        println!();
    }
    println!("];");
}
