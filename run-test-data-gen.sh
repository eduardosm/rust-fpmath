#!/bin/sh
set -e

DATA_DIR="$(dirname "$(readlink -f "$0")")/fpmath-tests/data"
mkdir -p "$DATA_DIR"

cargo run -p fpmath-tests --release -- --generate "$@"
