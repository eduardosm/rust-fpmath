#!/bin/sh
set -e

cargo run -p generator -- "$@"
