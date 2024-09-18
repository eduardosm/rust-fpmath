#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

target="x86_64-unknown-none"
features_array=("" "soft-float")

for features in "${features_array[@]}"; do
  begin_group "Build fpmath, target=\"$target\", features=\"$features\""
  cargo build --frozen -p fpmath --target "$target" --no-default-features --features "$features"
  end_group
done
