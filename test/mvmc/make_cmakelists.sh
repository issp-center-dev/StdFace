#!/bin/sh

for d in *
do
    if [ ! -d $d ]; then
        continue
    fi
    echo "add_test("
    echo "  NAME $d"
    echo "  COMMAND \${CMAKE_SOURCE_DIR}/test/mvmc/check.sh $d \${CMAKE_SOURCE_DIR}/test/mvmc"
    echo ")"
done
