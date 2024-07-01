#!/bin/sh

CMD=hwave_dry.out

for d in rpa* uhfk* uhfr*
do
  echo $d
  (
    cd $d
    if [ ! -f stan.in ]; then
      continue
    fi
    if [ -d ref ]; then
      rm -rf ref
    fi
    mkdir ref
    cd ref
    $CMD ../stan.in
  )
done
 
