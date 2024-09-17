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

begin_group "Generate test data"
./run-all-test-data-gen.sh
end_group

begin_group "Test x86_64"
cargo test --frozen --workspace "${features[@]}" --target x86_64-unknown-linux-gnu
end_group

begin_group "Test i586 debug"
cargo test --frozen --workspace "${features[@]}" --target i586-unknown-linux-gnu
end_group

begin_group "Test i586 release"
cargo test --frozen --workspace "${features[@]}" --target i586-unknown-linux-gnu --release
end_group
