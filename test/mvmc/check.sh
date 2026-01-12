#!/bin/sh
set -eu
test_item=${1:?}
base_dir=${2:?}
exec "${base_dir}/../common/check_case.sh" \
  "$test_item" "$base_dir" "StdFace.def" "../../../src/mvmc_dry.out"
