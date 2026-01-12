#!/bin/sh
set -eu
test_item=${1:?}
base_dir=${2:?}
exec "${base_dir}/../common/check_case.sh" \
  "$test_item" "$base_dir" "stan.in" "../../../src/hwave_dry.out"
