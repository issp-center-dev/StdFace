#!/bin/sh
set -eu
binary=${1:?}
workdir=${2:?}
mkdir -p "$workdir"

# Compare every on-site matrix element with h Sz + Gamma Sx + Gamma_y Sy,
# including the localized and itinerant sites, zero fields and general spin.
for s2 in 1 2; do
  for cells in 1 2; do
    for fields in zero z x y xyz; do
      h=0; gx=0; gy=0
      case "$fields" in
        z) h=0.2;;
        x) gx=0.4;;
        y) gy=0.6;;
        xyz) h=0.2; gx=0.4; gy=0.6;;
      esac
      rundir="$workdir/s${s2}-n${cells}-${fields}"
      mkdir -p "$rundir"
      (
        cd "$rundir"
        cat > stan.in <<EOF
model = "KondoGC"
method = "CG"
lattice = "fcc"
L = $cells
W = 1
Height = 1
2S = $s2
t = 0
J = 1
mu = 0
h = $h
Gamma = $gx
Gamma_y = $gy
EOF
        if ! "$binary" stan.in > run.log 2>&1; then
          cat run.log
          exit 1
        fi
        awk -v n="$cells" -v s2="$s2" -v h="$h" -v gx="$gx" -v gy="$gy" '
          function abs(x) { return x < 0 ? -x : x }
          NR > 5 {
            maxspin = $1 < n ? s2 : 1
            if (NF != 6 || $1 != $3 || $1 < 0 || $1 >= 2*n ||
                $2 < 0 || $2 > maxspin || $4 < 0 || $4 > maxspin) bad=1
            re[$1,$2,$4] += $5
            im[$1,$2,$4] += $6
          }
          END {
            for (site=0; site<2*n; site++) {
              twiceS = site < n ? s2 : 1
              for (a=0; a<=twiceS; a++) for (b=0; b<=twiceS; b++) {
                er=0; ei=0
                if (a==b) er=h*(twiceS/2-a)
                if (abs(a-b)==1) {
                  k = a>b ? a : b
                  amplitude=sqrt(k*(twiceS-k+1))/2
                  er=gx*amplitude
                  ei=(a>b ? gy : -gy)*amplitude
                }
                if (abs(re[site,a,b]-er)>1e-10 || abs(im[site,a,b]-ei)>1e-10) {
                  printf "site=%d spin=(%d,%d): got (%g,%g), expected (%g,%g)\n", site,a,b,re[site,a,b],im[site,a,b],er,ei
                  bad=1
                }
              }
            }
            exit bad
          }
        ' trans.def
      )
    done
  done
done
