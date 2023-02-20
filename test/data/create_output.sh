#!/bin/bash
cd "$(dirname "$0")"
set -Eeuo pipefail

if ! which valik 2> /dev/null; then
    echo "valik not found. please compile it and add it the path to it to the PATH variable"
    exit 255
fi

export VALIK_TMP=tmp/valik/my_dir
export VALIK_STELLAR=echo

./split/cli_test_output.sh

./build/cli_test_output.sh

./search/cli_test_output.sh

./consolidate/cli_test_output.sh
