#!/bin/sh
set -eu
test_item=${1:?}
base_dir=${2:?}
mode=${3:-"base"}

if [ "$mode" = "base" ]; then
    exec "${base_dir}/../common/check_case.sh" \
	 "$test_item" "$base_dir" "StdFace.def" "../../../src/mvmc_dry.out"
elif [ "$mode" = "python" ]; then
    exec "${base_dir}/../common/check_case.sh" \
	 "$test_item" "$base_dir" "StdFace.def" "python3 -mstdface --solver mVMC"
else
    /usr/bin/false
fi
