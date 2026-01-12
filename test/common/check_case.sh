#!/bin/sh
# Common CTest helper for StdFace:
# - Copy case directory into working directory
# - Run dry-run generator
# - Compare outputs against ref/ using diff or almost_diff.py
#
# Usage:
#   check_case.sh <test_item> <base_dir> <input_file> <default_cmd>
#
# Override command:
#   DRY_CMD="python3 -m stdface.hphi_dry" ctest -V -R <testname>

set -eu

test_item=${1:?}
base_dir=${2:?}
input_file=${3:?}
default_cmd=${4:?}

# Allow caller override (may include args)
DRY_CMD="${DRY_CMD:-$default_cmd}"

almost_diff_py="${base_dir}/../almost_diff.py"

# Avoid collisions across repeated runs
if [ -d "${test_item}.bak" ]; then
    rm -rf "${test_item}.bak"
fi
if [ -d "${test_item}" ]; then
    mv "${test_item}" "${test_item}.bak"
fi

# Copy case directory into working directory (CTest sets WORKING_DIRECTORY)
cp -rp "${base_dir}/${test_item}" .

cd "${test_item}"

# Run generator (log to run.log to match legacy behavior)
sh -c "$DRY_CMD \"$input_file\"" > run.log 2>&1

# Compare all files listed under ref/
# Keep legacy semantics (ls-based list)
file_list="$(cd ref && ls)"

err=""

for f in ${file_list}; do
    {
        echo "testing ${f}..."
        if diff "ref/$f" "$f" \
          || python3 "$almost_diff_py" "ref/$f" "$f"; then
            echo "ok"
        else
            echo "$f: reference and result differ"
            err="${err} ${f}"
        fi
    } >> run.log 2>&1
done

if [ -n "${err}" ]; then
    echo "error found in:${err}"
    exit 1
fi

exit 0
