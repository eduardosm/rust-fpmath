pub(super) fn expect_0_args(args: &[&str]) -> Result<(), String> {
    if args.is_empty() {
        Ok(())
    } else {
        Err(format!("expected 0 arguments, found {}", args.len()))
    }
}

pub(super) fn parse_1_arg<T0>(args: &[&str]) -> Result<T0, String>
where
    // Use associated type bounds (FromStr<Err: Display>)
    // when MSRV is 1.79 or higher
    T0: std::str::FromStr,
    <T0 as std::str::FromStr>::Err: std::fmt::Display,
{
    if let [arg1] = args {
        let v1 = arg1
            .parse()
            .map_err(|e| format!("failed to parse first argument {arg1:?}: {e}"))?;
        Ok(v1)
    } else {
        Err(format!("expected 1 argument, found {}", args.len()))
    }
}

pub(super) fn parse_2_args<T0, T1>(args: &[&str]) -> Result<(T0, T1), String>
where
    T0: std::str::FromStr,
    <T0 as std::str::FromStr>::Err: std::fmt::Display,
    T1: std::str::FromStr,
    <T1 as std::str::FromStr>::Err: std::fmt::Display,
{
    if let [arg1, arg2] = args {
        let v1 = arg1
            .parse()
            .map_err(|e| format!("failed to parse first argument {arg1:?}: {e}"))?;
        let v2 = arg2
            .parse()
            .map_err(|e| format!("failed to parse second argument {arg2:?}: {e}"))?;
        Ok((v1, v2))
    } else {
        Err(format!("expected 2 arguments, found {}", args.len()))
    }
}

pub(super) fn parse_3_args<T0, T1, T2>(args: &[&str]) -> Result<(T0, T1, T2), String>
where
    T0: std::str::FromStr,
    <T0 as std::str::FromStr>::Err: std::fmt::Display,
    T1: std::str::FromStr,
    <T1 as std::str::FromStr>::Err: std::fmt::Display,
    T2: std::str::FromStr,
    <T2 as std::str::FromStr>::Err: std::fmt::Display,
{
    if let [arg1, arg2, arg3] = args {
        let v1 = arg1
            .parse()
            .map_err(|e| format!("failed to parse first argument {arg1:?}: {e}"))?;
        let v2 = arg2
            .parse()
            .map_err(|e| format!("failed to parse second argument {arg2:?}: {e}"))?;
        let v3 = arg3
            .parse()
            .map_err(|e| format!("failed to parse third argument {arg3:?}: {e}"))?;
        Ok((v1, v2, v3))
    } else {
        Err(format!("expected 3 arguments, found {}", args.len()))
    }
}

pub(super) fn parse_4_args<T0, T1, T2, T3>(args: &[&str]) -> Result<(T0, T1, T2, T3), String>
where
    T0: std::str::FromStr,
    <T0 as std::str::FromStr>::Err: std::fmt::Display,
    T1: std::str::FromStr,
    <T1 as std::str::FromStr>::Err: std::fmt::Display,
    T2: std::str::FromStr,
    <T2 as std::str::FromStr>::Err: std::fmt::Display,
    T3: std::str::FromStr,
    <T3 as std::str::FromStr>::Err: std::fmt::Display,
{
    if let [arg1, arg2, arg3, arg4] = args {
        let v1 = arg1
            .parse()
            .map_err(|e| format!("failed to parse first argument {arg1:?}: {e}"))?;
        let v2 = arg2
            .parse()
            .map_err(|e| format!("failed to parse second argument {arg2:?}: {e}"))?;
        let v3 = arg3
            .parse()
            .map_err(|e| format!("failed to parse third argument {arg3:?}: {e}"))?;
        let v4 = arg4
            .parse()
            .map_err(|e| format!("failed to parse fourth argument {arg4:?}: {e}"))?;
        Ok((v1, v2, v3, v4))
    } else {
        Err(format!("expected 4 arguments, found {}", args.len()))
    }
}
