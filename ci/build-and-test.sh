#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

export RUSTDOCFLAGS="-D warnings"

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Build"
cargo build --frozen --workspace --all-targets
end_group

begin_group "Doc"
cargo doc --frozen --workspace
end_group

begin_group "Test"
cargo test --frozen --workspace
end_group
