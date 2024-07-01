#!/bin/sh

CMD=mvmc_dry.out
#CMD=vmcdry.out
INFILE=StdFace.def

for d in *
do
  if [ ! -d $d ]; then
    continue
  fi
  echo $d
  (
    cd $d
    if [ ! -f $INFILE ]; then
      continue
    fi
    if [ -d ref ]; then
      rm -rf ref
    fi
    mkdir ref
    cd ref
    $CMD ../$INFILE
  )
done
 
