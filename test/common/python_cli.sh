#!/bin/sh
# Stand-in for a C dry-run binary in python mode: runs the Python CLI for the
# solver named in STDFACE_SOLVER (HPhi, mVMC, UHF, HWAVE).
# The C generators exit with 255 on an input error and the scripts written
# for them check for that; the Python CLI exits with 1. An uncaught Python
# exception is not an input error and keeps exit status 1.
out=$(stdface --solver "${STDFACE_SOLVER:?}" "$@" 2>&1)
status=$?
printf '%s\n' "$out"
[ $status -eq 0 ] && exit 0
case "$out" in
    *"Traceback (most recent call last)"*) exit 1 ;;
esac
exit 255
