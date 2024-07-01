#!/bin/sh

CMD=uhf_dry.out
INFILE=stan.in

for d in *
do
  if [ ! -d $d ]; then
    continue
  fi
  if [ ! -f $d/$INFILE ]; then
    continue
  fi 
  echo $d
  (
    cd $d
    if [ -d ref ]; then
      rm -rf ref
    fi
    mkdir ref
    cd ref
    $CMD ../$INFILE
  )
done
 
