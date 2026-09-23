#!/usr/bin/env bash
set -euo pipefail

. ci/utils.sh

begin_group "Install Julia packages"
julia -e 'import Pkg; Pkg.add(name="Remez", version="0.1.1"); Pkg.add("SpecialFunctions")'
end_group

begin_group "Fetch dependencies"
cargo fetch --locked
end_group

begin_group "Run generator"
cargo run --frozen -p generator -- src
end_group

begin_group "Show diff"
git diff --exit-code
end_group
