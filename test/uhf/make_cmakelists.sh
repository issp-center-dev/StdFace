#!/bin/sh

INFILE=stan.in

for d in *
do
    if [ ! -d $d ]; then
        continue
    fi
    if [ ! -f $d/$INFILE ]; then
        continue
    fi
    echo "add_test("
    echo "  NAME $d"
    echo "  COMMAND \${CMAKE_SOURCE_DIR}/test/uhf/check.sh $d \${CMAKE_SOURCE_DIR}/test/hwave"
    echo ")"
done
