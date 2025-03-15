#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Install Julia packages"
julia -e 'import Pkg; Pkg.add("Remez"); Pkg.add("SpecialFunctions")'
end_group

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Run generator"
cargo run --frozen -p generator -- src
end_group

begin_group "Build after generation"
cargo build --frozen -p fpmath
end_group

begin_group "Show diff"
git diff
end_group
