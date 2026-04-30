#!/bin/sh

test_item=$1
base_dir=$2
mode=${3:-"base"}

if [ "$mode" = "base" ]; then
    DRY_BIN=${DRY_BIN:-../../../src/hphi_dry.out}
elif [ "$mode" = "python" ]; then
    DRY_BIN=${DRY_BIN:-"python3 -mstdface --solver HPhi"}
else
    DRY_BIN="/usr/bin/false"
fi

if [ -d ${test_item} ]; then
    mv ${test_item} ${test_item}.bak
fi

cp -rp ${base_dir}/${test_item} .
cd ${test_item}

${DRY_BIN} stan.in > run.log 2>&1
status=$?

if [ $status -eq 0 ]; then
    echo "expected failure but exit code was 0"
    exit 1
fi

exit 0
