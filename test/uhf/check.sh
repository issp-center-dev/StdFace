#!/bin/sh

test_item=$1
base_dir=$2

DRY_BIN=../../../src/uhf_dry.out

almost_diff="python3 ${base_dir}/../almost_diff.py"

if [ -d ${test_item} ]; then
    mv ${test_item} ${test_item}.bak
fi

cp -rp ${base_dir}/${test_item} .

cd ${test_item}

${DRY_BIN} stan.in > run.log 2>&1

file_list=`(cd ref && ls)`

err=""

for f in ${file_list}
do
    echo "testing ${f}..."
    if diff ref/$f $f || $almost_diff ref/$f $f ; then
	echo "ok"
    else
	echo "$f: reference and result differ"
	err="$err $f"
    fi
done >> run.log 2>&1

if [ -n "$err" ]; then
    echo "error found in:$err"
    exit 1
fi

exit 0
