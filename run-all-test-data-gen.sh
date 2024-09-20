#!/bin/sh
set -e

DATA_DIR="$(dirname "$(readlink -f "$0")")/fpmath-tests/data"
if [ -e "$DATA_DIR" ]; then
    rm -r "$DATA_DIR"
fi

./run-test-data-gen.sh --all
