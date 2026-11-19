#!/bin/sh
set -eu

test_item=${1:?}
base_dir=${2:?}
dry_cmd=${3:?}

if [ -d "${test_item}.bak" ]; then
    rm -rf "${test_item}.bak"
fi
if [ -d "${test_item}" ]; then
    mv "${test_item}" "${test_item}.bak"
fi

cp -rp "${base_dir}/${test_item}" .
cd "${test_item}"

"${dry_cmd}" StdFace.def > run.log 2>&1

if ! grep -q "^NSROptItrStep  2$" modpara.def; then
    echo "NSROptItrStep was not set to 2"
    exit 1
fi

if ! grep -q "^NSROptItrSmp   1$" modpara.def; then
    echo "NSROptItrSmp default was not clamped to 1"
    exit 1
fi
