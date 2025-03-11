pub(crate) fn gen_frac_2_pi_large() {
    let num_words = 66;
    let mut tmp = 2u8 / rug::Float::with_val(num_words * 24, rug::float::Constant::Pi);

    println!("const FRAC_2_PI_LARGE: &[u32] = &[");
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            print!("    ");
        } else {
            print!(" ");
        }
        tmp <<= 24;
        let word = tmp
            .to_u32_saturating_round(rug::float::Round::Zero)
            .unwrap();
        print!("0x{word:06X},");
        tmp -= word;
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
    let mut tmp = rug::Float::with_val(num_words * 24, rug::float::Constant::Pi) / 2u8;

    println!("const FRAC_PI_2_MEDIUM: &[u32] = &[");
    tmp /= 2u8;
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            print!("    ");
        } else {
            print!(" ");
        }
        tmp <<= 24;
        let word = tmp
            .to_u32_saturating_round(rug::float::Round::Zero)
            .unwrap();
        print!("0x{word:06X},");
        tmp -= word;
        if (i % words_per_line) == (words_per_line - 1) {
            println!();
        }
    }
    if (num_words % words_per_line) != 0 {
        println!();
    }
    println!("];");
}
