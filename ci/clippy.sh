#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Run clippy workspace"
cargo clippy --frozen --workspace --all-targets  -- -D warnings
end_group

features_array=("" "soft-float")

for features in "${features_array[@]}"; do
  begin_group "Run clippy fpmath, features=\"$features\""
  cargo clippy --frozen -p fpmath --all-targets --no-default-features --features "$features" -- -D warnings
  end_group
done
