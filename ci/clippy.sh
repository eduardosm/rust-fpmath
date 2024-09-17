#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Run clippy workspace"
cargo clippy --frozen --workspace --all-targets  -- -D warnings
end_group

begin_group "Run clippy fpmath without soft-float"
cargo clippy --frozen -p fpmath --all-targets -- -D warnings
end_group

begin_group "Run clippy fpmath with soft-float"
cargo clippy --frozen -p fpmath --all-targets --features "soft-float" -- -D warnings
end_group
