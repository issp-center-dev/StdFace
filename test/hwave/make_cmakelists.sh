#!/bin/sh

for d in uhfr* uhfk* rpa*
do
    echo "add_test("
    echo "  NAME $d"
    echo "  COMMAND \${CMAKE_SOURCE_DIR}/test/hwave/check.sh $d \${CMAKE_SOURCE_DIR}/test/hwave"
    echo ")"
done
