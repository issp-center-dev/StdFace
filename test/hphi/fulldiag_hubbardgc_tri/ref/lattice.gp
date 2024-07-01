#set terminal pdf color enhanced \
#dashed dl 1.0 size 20.0cm, 20.0cm 
#set output "lattice.pdf"
set xrange [-2.000000: 4.598076]
set yrange [-2.000000: 4.598076]
set size square
unset key
unset tics
unset border
set style line 1 lc 1 lt 1
set style line 2 lc 5 lt 1
set style line 3 lc 0 lt 1
set arrow from 0.000000, 0.000000 to 1.500000, 0.866025 nohead front ls 3
set arrow from 1.500000, 0.866025 to 1.500000, 2.598076 nohead front ls 3
set arrow from 1.500000, 2.598076 to 0.000000, 1.732051 nohead front ls 3
set arrow from 0.000000, 1.732051 to 0.000000, 0.000000 nohead front ls 3
set label "0" at 0.000000, 0.000000 center front
set label "1" at -1.000000, 0.000000 center front
set arrow from 0.000000, 0.000000 to -1.000000, 0.000000 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "2" at 1.000000, 0.000000 center front
set arrow from 0.000000, 0.000000 to 1.000000, 0.000000 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "2" at -0.500000, -0.866025 center front
set arrow from 0.000000, 0.000000 to -0.500000, -0.866025 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "1" at 0.500000, 0.866025 center front
set arrow from 0.000000, 0.000000 to 0.500000, 0.866025 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "2" at -0.500000, 0.866025 center front
set arrow from 0.000000, 0.000000 to -0.500000, 0.866025 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "1" at 0.500000, -0.866025 center front
set arrow from 0.000000, 0.000000 to 0.500000, -0.866025 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "0" at -1.500000, 0.866025 center front
set arrow from 0.000000, 0.000000 to -1.500000, 0.866025 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "0" at 1.500000, -0.866025 center front
set arrow from 0.000000, 0.000000 to 1.500000, -0.866025 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "0" at -1.500000, -0.866025 center front
set arrow from 0.000000, 0.000000 to -1.500000, -0.866025 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "0" at 1.500000, 0.866025 center front
set arrow from 0.000000, 0.000000 to 1.500000, 0.866025 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "0" at 0.000000, -1.732051 center front
set arrow from 0.000000, 0.000000 to 0.000000, -1.732051 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "0" at 0.000000, 1.732051 center front
set arrow from 0.000000, 0.000000 to 0.000000, 1.732051 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "2" at -2.000000, 0.000000 center front
set label "0" at 0.000000, 0.000000 center front
set label "1" at 2.000000, 0.000000 center front
set label "0" at 0.000000, 0.000000 center front
set label "1" at -1.000000, -1.732051 center front
set label "0" at 0.000000, 0.000000 center front
set label "2" at 1.000000, 1.732051 center front
set label "0" at 0.000000, 0.000000 center front
set label "1" at -1.000000, 1.732051 center front
set label "0" at 0.000000, 0.000000 center front
set label "2" at 1.000000, -1.732051 center front
set label "1" at 0.500000, 0.866025 center front
set label "2" at -0.500000, 0.866025 center front
set arrow from 0.500000, 0.866025 to -0.500000, 0.866025 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "0" at 1.500000, 0.866025 center front
set arrow from 0.500000, 0.866025 to 1.500000, 0.866025 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "0" at 0.000000, 0.000000 center front
set arrow from 0.500000, 0.866025 to 0.000000, 0.000000 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "2" at 1.000000, 1.732051 center front
set arrow from 0.500000, 0.866025 to 1.000000, 1.732051 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "0" at 0.000000, 1.732051 center front
set arrow from 0.500000, 0.866025 to 0.000000, 1.732051 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "2" at 1.000000, 0.000000 center front
set arrow from 0.500000, 0.866025 to 1.000000, 0.000000 nohead ls 1
set label "1" at 0.500000, 0.866025 center front
set label "1" at -1.000000, 1.732051 center front
set arrow from 0.500000, 0.866025 to -1.000000, 1.732051 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "1" at 2.000000, 0.000000 center front
set arrow from 0.500000, 0.866025 to 2.000000, 0.000000 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "1" at -1.000000, 0.000000 center front
set arrow from 0.500000, 0.866025 to -1.000000, 0.000000 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "1" at 2.000000, 1.732051 center front
set arrow from 0.500000, 0.866025 to 2.000000, 1.732051 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "1" at 0.500000, -0.866025 center front
set arrow from 0.500000, 0.866025 to 0.500000, -0.866025 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "1" at 0.500000, 2.598076 center front
set arrow from 0.500000, 0.866025 to 0.500000, 2.598076 nohead ls 2
set label "1" at 0.500000, 0.866025 center front
set label "0" at -1.500000, 0.866025 center front
set label "1" at 0.500000, 0.866025 center front
set label "2" at 2.500000, 0.866025 center front
set label "1" at 0.500000, 0.866025 center front
set label "2" at -0.500000, -0.866025 center front
set label "1" at 0.500000, 0.866025 center front
set label "0" at 1.500000, 2.598076 center front
set label "1" at 0.500000, 0.866025 center front
set label "2" at -0.500000, 2.598076 center front
set label "1" at 0.500000, 0.866025 center front
set label "0" at 1.500000, -0.866025 center front
set label "2" at 1.000000, 1.732051 center front
set label "0" at 0.000000, 1.732051 center front
set arrow from 1.000000, 1.732051 to 0.000000, 1.732051 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "1" at 2.000000, 1.732051 center front
set arrow from 1.000000, 1.732051 to 2.000000, 1.732051 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "1" at 0.500000, 0.866025 center front
set arrow from 1.000000, 1.732051 to 0.500000, 0.866025 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "0" at 1.500000, 2.598076 center front
set arrow from 1.000000, 1.732051 to 1.500000, 2.598076 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "1" at 0.500000, 2.598076 center front
set arrow from 1.000000, 1.732051 to 0.500000, 2.598076 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "0" at 1.500000, 0.866025 center front
set arrow from 1.000000, 1.732051 to 1.500000, 0.866025 nohead ls 1
set label "2" at 1.000000, 1.732051 center front
set label "2" at -0.500000, 2.598076 center front
set arrow from 1.000000, 1.732051 to -0.500000, 2.598076 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "2" at 2.500000, 0.866025 center front
set arrow from 1.000000, 1.732051 to 2.500000, 0.866025 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "2" at -0.500000, 0.866025 center front
set arrow from 1.000000, 1.732051 to -0.500000, 0.866025 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "2" at 2.500000, 2.598076 center front
set arrow from 1.000000, 1.732051 to 2.500000, 2.598076 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "2" at 1.000000, 0.000000 center front
set arrow from 1.000000, 1.732051 to 1.000000, 0.000000 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "2" at 1.000000, 3.464102 center front
set arrow from 1.000000, 1.732051 to 1.000000, 3.464102 nohead ls 2
set label "2" at 1.000000, 1.732051 center front
set label "1" at -1.000000, 1.732051 center front
set label "2" at 1.000000, 1.732051 center front
set label "0" at 3.000000, 1.732051 center front
set label "2" at 1.000000, 1.732051 center front
set label "0" at 0.000000, 0.000000 center front
set label "2" at 1.000000, 1.732051 center front
set label "1" at 2.000000, 3.464102 center front
set label "2" at 1.000000, 1.732051 center front
set label "0" at 0.000000, 3.464102 center front
set label "2" at 1.000000, 1.732051 center front
set label "1" at 2.000000, 0.000000 center front
plot '-' w d lc 7
0.0 0.0
end
pause -1
