#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Test i586 debug"
cargo test --frozen --workspace --target i586-unknown-linux-gnu
end_group

begin_group "Test i586 release"
cargo test --frozen --workspace --target i586-unknown-linux-gnu --release
end_group
