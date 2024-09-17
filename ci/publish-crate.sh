#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

export CARGO_REGISTRY_TOKEN="$CRATES_IO_TOKEN"

crate="fpmath"

begin_group "Publish $crate"
cargo publish --locked -p "$crate" --no-verify
end_group
