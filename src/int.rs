use crate::traits;

macro_rules! impl_cast_from {
    ($s:ty as $d:ty) => {
        impl traits::CastFrom<$s> for $d {
            #[allow(trivial_numeric_casts)]
            #[inline]
            fn cast_from(value: $s) -> $d {
                value as $d
            }
        }
    };
}

macro_rules! impl_int {
    ($t:ty) => {
        impl_cast_from!($t as u8);
        impl_cast_from!($t as i8);
        impl_cast_from!($t as u16);
        impl_cast_from!($t as i16);
        impl_cast_from!($t as u32);
        impl_cast_from!($t as i32);
        impl_cast_from!($t as u64);
        impl_cast_from!($t as i64);
        impl_cast_from!($t as f32);
        impl_cast_from!($t as f64);

        impl traits::Int for $t {
            const ZERO: Self = 0;
            const ONE: Self = 1;
            const TWO: Self = 2;

            const MAX: Self = <$t>::MAX;
        }
    };
}

macro_rules! impl_uint {
    ($t:ty) => {
        impl_int!($t);
        impl traits::UInt for $t {}
    };
}

macro_rules! impl_sint {
    ($t:ty) => {
        impl_int!($t);
        impl traits::SInt for $t {}
    };
}

impl_uint!(u8);
impl_uint!(u16);
impl_uint!(u32);
impl_uint!(u64);

impl_sint!(i8);
impl_sint!(i16);
impl_sint!(i32);
impl_sint!(i64);
