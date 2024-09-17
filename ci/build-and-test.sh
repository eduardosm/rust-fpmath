#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

export RUSTDOCFLAGS="-D warnings"

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Build"
cargo build --workspace --all-targets --frozen
end_group

begin_group "Doc"
cargo doc --workspace --frozen
end_group

begin_group "Generate test data"
./run-all-test-data-gen.sh
end_group

begin_group "Test x86_64"
cargo test --workspace --target x86_64-unknown-linux-gnu --frozen
end_group

begin_group "Test i586 debug"
cargo test --workspace --target i586-unknown-linux-gnu --frozen
end_group

begin_group "Test i586 release"
cargo test --workspace --target i586-unknown-linux-gnu --release --frozen
end_group
