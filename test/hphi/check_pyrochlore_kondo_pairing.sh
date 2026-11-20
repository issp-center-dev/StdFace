#!/bin/sh
set -eu
binary=${1:?}
workdir=${2:?}
mkdir -p "$workdir"

# The Kondo coupling must pair each itinerant site with the localized spin of
# its own sublattice: localized site s couples to itinerant site s + n only,
# and every localized site is coupled exactly once.
for cells in 1 2; do
  rundir="$workdir/n${cells}"
  mkdir -p "$rundir"
  (
    cd "$rundir"
    cat > stan.in <<EOF
model = "KondoGC"
method = "CG"
lattice = "pyrochlore"
L = $cells
W = 1
Height = 1
t = 0
J = 1
EOF
    if ! "$binary" stan.in > run.log 2>&1; then
      cat run.log
      exit 1
    fi
    for file in exchange.def hund.def coulombinter.def; do
      awk -v n="$((4 * cells))" -v file="$file" '
        NR > 5 {
          lo = $1 < $2 ? $1 : $2
          hi = $1 < $2 ? $2 : $1
          if (NF != 3 || lo < 0 || lo >= n || hi != lo + n) {
            printf "%s: unexpected pair (%d,%d)\n", file, $1, $2
            bad=1
          }
          seen[lo]++
        }
        END {
          for (site=0; site<n; site++) {
            if (seen[site] != 1) {
              printf "%s: localized site %d is coupled %d times\n", file, site, seen[site]
              bad=1
            }
          }
          exit bad
        }
      ' "$file"
    done
  )
done
