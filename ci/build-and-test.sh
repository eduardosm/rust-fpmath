#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

export RUSTDOCFLAGS="-D warnings"

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

features=(--features "soft-float")

begin_group "Build"
cargo build --frozen --workspace "${features[@]}" --all-targets
end_group

begin_group "Doc"
cargo doc --frozen --workspace "${features[@]}"
end_group

begin_group "Test"
cargo test --frozen --workspace "${features[@]}"
end_group
