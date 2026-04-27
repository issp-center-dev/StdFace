#set terminal pdf color enhanced \
#dashed dl 1.0 size 20.0cm, 20.0cm 
#set output "lattice.pdf"
set xrange [-2.000000: 5.000000]
set yrange [-2.000000: 5.000000]
set size square
unset key
unset tics
unset border
set style line 1 lc 1 lt 1
set style line 2 lc 5 lt 1
set style line 3 lc 0 lt 1
set arrow from 0.000000, 0.000000 to 2.000000, 0.000000 nohead front ls 3
set arrow from 2.000000, 0.000000 to 3.000000, 1.732051 nohead front ls 3
set arrow from 3.000000, 1.732051 to 1.000000, 1.732051 nohead front ls 3
set arrow from 1.000000, 1.732051 to 0.000000, 0.000000 nohead front ls 3
set label "13" at 0.500000, 0.000000 center front
set label "12" at 0.000000, 0.000000 center front
set arrow from 0.500000, 0.000000 to 0.000000, 0.000000 nohead ls 1
set label "12" at 0.000000, 0.000000 center front
set label "13" at 0.500000, 0.000000 center front
set arrow from 0.000000, 0.000000 to 0.500000, 0.000000 nohead ls 1
set label "14" at 0.250000, 0.433013 center front
set label "12" at 0.000000, 0.000000 center front
set arrow from 0.250000, 0.433013 to 0.000000, 0.000000 nohead ls 1
set label "12" at 0.000000, 0.000000 center front
set label "14" at 0.250000, 0.433013 center front
set arrow from 0.000000, 0.000000 to 0.250000, 0.433013 nohead ls 1
set label "14" at 0.250000, 0.433013 center front
set label "13" at 0.500000, 0.000000 center front
set arrow from 0.250000, 0.433013 to 0.500000, 0.000000 nohead ls 1
set label "13" at 0.500000, 0.000000 center front
set label "14" at 0.250000, 0.433013 center front
set arrow from 0.500000, 0.000000 to 0.250000, 0.433013 nohead ls 1
set label "12" at 0.000000, 0.000000 center front
set label "16" at -0.500000, 0.000000 center front
set arrow from 0.000000, 0.000000 to -0.500000, 0.000000 nohead ls 1
set label "13" at 0.500000, 0.000000 center front
set label "15" at 1.000000, 0.000000 center front
set arrow from 0.500000, 0.000000 to 1.000000, 0.000000 nohead ls 1
set label "12" at 0.000000, 0.000000 center front
set label "20" at -0.250000, -0.433013 center front
set arrow from 0.000000, 0.000000 to -0.250000, -0.433013 nohead ls 1
set label "14" at 0.250000, 0.433013 center front
set label "18" at 0.500000, 0.866025 center front
set arrow from 0.250000, 0.433013 to 0.500000, 0.866025 nohead ls 1
set label "14" at 0.250000, 0.433013 center front
set label "22" at 0.000000, 0.866025 center front
set arrow from 0.250000, 0.433013 to 0.000000, 0.866025 nohead ls 1
set label "13" at 0.500000, 0.000000 center front
set label "23" at 0.750000, -0.433013 center front
set arrow from 0.500000, 0.000000 to 0.750000, -0.433013 nohead ls 1
set label "12" at 0.000000, 0.000000 center front
set label "17" at -0.750000, 0.433013 center front
set arrow from 0.000000, 0.000000 to -0.750000, 0.433013 nohead ls 2
set label "14" at 0.250000, 0.433013 center front
set label "15" at 1.000000, 0.000000 center front
set arrow from 0.250000, 0.433013 to 1.000000, 0.000000 nohead ls 2
set label "14" at 0.250000, 0.433013 center front
set label "16" at -0.500000, 0.000000 center front
set arrow from 0.250000, 0.433013 to -0.500000, 0.000000 nohead ls 2
set label "13" at 0.500000, 0.000000 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 0.500000, 0.000000 to 1.250000, 0.433013 nohead ls 2
set label "12" at 0.000000, 0.000000 center front
set label "19" at 0.000000, -0.866025 center front
set arrow from 0.000000, 0.000000 to 0.000000, -0.866025 nohead ls 2
set label "13" at 0.500000, 0.000000 center front
set label "18" at 0.500000, 0.866025 center front
set arrow from 0.500000, 0.000000 to 0.500000, 0.866025 nohead ls 2
set label "13" at 0.500000, 0.000000 center front
set label "20" at -0.250000, -0.433013 center front
set arrow from 0.500000, 0.000000 to -0.250000, -0.433013 nohead ls 2
set label "14" at 0.250000, 0.433013 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 0.250000, 0.433013 to 1.000000, 0.866025 nohead ls 2
set label "14" at 0.250000, 0.433013 center front
set label "21" at -0.500000, 0.866025 center front
set arrow from 0.250000, 0.433013 to -0.500000, 0.866025 nohead ls 2
set label "12" at 0.000000, 0.000000 center front
set label "23" at 0.750000, -0.433013 center front
set arrow from 0.000000, 0.000000 to 0.750000, -0.433013 nohead ls 2
set label "13" at 0.500000, 0.000000 center front
set label "21" at 0.500000, -0.866025 center front
set arrow from 0.500000, 0.000000 to 0.500000, -0.866025 nohead ls 2
set label "12" at 0.000000, 0.000000 center front
set label "22" at 0.000000, 0.866025 center front
set arrow from 0.000000, 0.000000 to 0.000000, 0.866025 nohead ls 2
set label "16" at 1.500000, 0.000000 center front
set label "15" at 1.000000, 0.000000 center front
set arrow from 1.500000, 0.000000 to 1.000000, 0.000000 nohead ls 1
set label "15" at 1.000000, 0.000000 center front
set label "16" at 1.500000, 0.000000 center front
set arrow from 1.000000, 0.000000 to 1.500000, 0.000000 nohead ls 1
set label "17" at 1.250000, 0.433013 center front
set label "15" at 1.000000, 0.000000 center front
set arrow from 1.250000, 0.433013 to 1.000000, 0.000000 nohead ls 1
set label "15" at 1.000000, 0.000000 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 1.000000, 0.000000 to 1.250000, 0.433013 nohead ls 1
set label "17" at 1.250000, 0.433013 center front
set label "16" at 1.500000, 0.000000 center front
set arrow from 1.250000, 0.433013 to 1.500000, 0.000000 nohead ls 1
set label "16" at 1.500000, 0.000000 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 1.500000, 0.000000 to 1.250000, 0.433013 nohead ls 1
set label "15" at 1.000000, 0.000000 center front
set label "13" at 0.500000, 0.000000 center front
set arrow from 1.000000, 0.000000 to 0.500000, 0.000000 nohead ls 1
set label "16" at 1.500000, 0.000000 center front
set label "12" at 2.000000, 0.000000 center front
set arrow from 1.500000, 0.000000 to 2.000000, 0.000000 nohead ls 1
set label "15" at 1.000000, 0.000000 center front
set label "23" at 0.750000, -0.433013 center front
set arrow from 1.000000, 0.000000 to 0.750000, -0.433013 nohead ls 1
set label "17" at 1.250000, 0.433013 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 1.250000, 0.433013 to 1.500000, 0.866025 nohead ls 1
set label "17" at 1.250000, 0.433013 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 1.250000, 0.433013 to 1.000000, 0.866025 nohead ls 1
set label "16" at 1.500000, 0.000000 center front
set label "20" at 1.750000, -0.433013 center front
set arrow from 1.500000, 0.000000 to 1.750000, -0.433013 nohead ls 1
set label "15" at 1.000000, 0.000000 center front
set label "14" at 0.250000, 0.433013 center front
set arrow from 1.000000, 0.000000 to 0.250000, 0.433013 nohead ls 2
set label "17" at 1.250000, 0.433013 center front
set label "12" at 2.000000, 0.000000 center front
set arrow from 1.250000, 0.433013 to 2.000000, 0.000000 nohead ls 2
set label "17" at 1.250000, 0.433013 center front
set label "13" at 0.500000, 0.000000 center front
set arrow from 1.250000, 0.433013 to 0.500000, 0.000000 nohead ls 2
set label "16" at 1.500000, 0.000000 center front
set label "14" at 2.250000, 0.433013 center front
set arrow from 1.500000, 0.000000 to 2.250000, 0.433013 nohead ls 2
set label "15" at 1.000000, 0.000000 center front
set label "22" at 1.000000, -0.866025 center front
set arrow from 1.000000, 0.000000 to 1.000000, -0.866025 nohead ls 2
set label "16" at 1.500000, 0.000000 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 1.500000, 0.000000 to 1.500000, 0.866025 nohead ls 2
set label "16" at 1.500000, 0.000000 center front
set label "23" at 0.750000, -0.433013 center front
set arrow from 1.500000, 0.000000 to 0.750000, -0.433013 nohead ls 2
set label "17" at 1.250000, 0.433013 center front
set label "22" at 2.000000, 0.866025 center front
set arrow from 1.250000, 0.433013 to 2.000000, 0.866025 nohead ls 2
set label "17" at 1.250000, 0.433013 center front
set label "18" at 0.500000, 0.866025 center front
set arrow from 1.250000, 0.433013 to 0.500000, 0.866025 nohead ls 2
set label "15" at 1.000000, 0.000000 center front
set label "20" at 1.750000, -0.433013 center front
set arrow from 1.000000, 0.000000 to 1.750000, -0.433013 nohead ls 2
set label "16" at 1.500000, 0.000000 center front
set label "18" at 1.500000, -0.866025 center front
set arrow from 1.500000, 0.000000 to 1.500000, -0.866025 nohead ls 2
set label "15" at 1.000000, 0.000000 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 1.000000, 0.000000 to 1.000000, 0.866025 nohead ls 2
set label "19" at 1.000000, 0.866025 center front
set label "18" at 0.500000, 0.866025 center front
set arrow from 1.000000, 0.866025 to 0.500000, 0.866025 nohead ls 1
set label "18" at 0.500000, 0.866025 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 0.500000, 0.866025 to 1.000000, 0.866025 nohead ls 1
set label "20" at 0.750000, 1.299038 center front
set label "18" at 0.500000, 0.866025 center front
set arrow from 0.750000, 1.299038 to 0.500000, 0.866025 nohead ls 1
set label "18" at 0.500000, 0.866025 center front
set label "20" at 0.750000, 1.299038 center front
set arrow from 0.500000, 0.866025 to 0.750000, 1.299038 nohead ls 1
set label "20" at 0.750000, 1.299038 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 0.750000, 1.299038 to 1.000000, 0.866025 nohead ls 1
set label "19" at 1.000000, 0.866025 center front
set label "20" at 0.750000, 1.299038 center front
set arrow from 1.000000, 0.866025 to 0.750000, 1.299038 nohead ls 1
set label "18" at 0.500000, 0.866025 center front
set label "22" at 0.000000, 0.866025 center front
set arrow from 0.500000, 0.866025 to 0.000000, 0.866025 nohead ls 1
set label "19" at 1.000000, 0.866025 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 1.000000, 0.866025 to 1.500000, 0.866025 nohead ls 1
set label "18" at 0.500000, 0.866025 center front
set label "14" at 0.250000, 0.433013 center front
set arrow from 0.500000, 0.866025 to 0.250000, 0.433013 nohead ls 1
set label "20" at 0.750000, 1.299038 center front
set label "12" at 1.000000, 1.732051 center front
set arrow from 0.750000, 1.299038 to 1.000000, 1.732051 nohead ls 1
set label "20" at 0.750000, 1.299038 center front
set label "16" at 0.500000, 1.732051 center front
set arrow from 0.750000, 1.299038 to 0.500000, 1.732051 nohead ls 1
set label "19" at 1.000000, 0.866025 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 1.000000, 0.866025 to 1.250000, 0.433013 nohead ls 1
set label "18" at 0.500000, 0.866025 center front
set label "23" at -0.250000, 1.299038 center front
set arrow from 0.500000, 0.866025 to -0.250000, 1.299038 nohead ls 2
set label "20" at 0.750000, 1.299038 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 0.750000, 1.299038 to 1.500000, 0.866025 nohead ls 2
set label "20" at 0.750000, 1.299038 center front
set label "22" at 0.000000, 0.866025 center front
set arrow from 0.750000, 1.299038 to 0.000000, 0.866025 nohead ls 2
set label "19" at 1.000000, 0.866025 center front
set label "23" at 1.750000, 1.299038 center front
set arrow from 1.000000, 0.866025 to 1.750000, 1.299038 nohead ls 2
set label "18" at 0.500000, 0.866025 center front
set label "13" at 0.500000, 0.000000 center front
set arrow from 0.500000, 0.866025 to 0.500000, 0.000000 nohead ls 2
set label "19" at 1.000000, 0.866025 center front
set label "12" at 1.000000, 1.732051 center front
set arrow from 1.000000, 0.866025 to 1.000000, 1.732051 nohead ls 2
set label "19" at 1.000000, 0.866025 center front
set label "14" at 0.250000, 0.433013 center front
set arrow from 1.000000, 0.866025 to 0.250000, 0.433013 nohead ls 2
set label "20" at 0.750000, 1.299038 center front
set label "13" at 1.500000, 1.732051 center front
set arrow from 0.750000, 1.299038 to 1.500000, 1.732051 nohead ls 2
set label "20" at 0.750000, 1.299038 center front
set label "15" at 0.000000, 1.732051 center front
set arrow from 0.750000, 1.299038 to 0.000000, 1.732051 nohead ls 2
set label "18" at 0.500000, 0.866025 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 0.500000, 0.866025 to 1.250000, 0.433013 nohead ls 2
set label "19" at 1.000000, 0.866025 center front
set label "15" at 1.000000, 0.000000 center front
set arrow from 1.000000, 0.866025 to 1.000000, 0.000000 nohead ls 2
set label "18" at 0.500000, 0.866025 center front
set label "16" at 0.500000, 1.732051 center front
set arrow from 0.500000, 0.866025 to 0.500000, 1.732051 nohead ls 2
set label "22" at 2.000000, 0.866025 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 2.000000, 0.866025 to 1.500000, 0.866025 nohead ls 1
set label "21" at 1.500000, 0.866025 center front
set label "22" at 2.000000, 0.866025 center front
set arrow from 1.500000, 0.866025 to 2.000000, 0.866025 nohead ls 1
set label "23" at 1.750000, 1.299038 center front
set label "21" at 1.500000, 0.866025 center front
set arrow from 1.750000, 1.299038 to 1.500000, 0.866025 nohead ls 1
set label "21" at 1.500000, 0.866025 center front
set label "23" at 1.750000, 1.299038 center front
set arrow from 1.500000, 0.866025 to 1.750000, 1.299038 nohead ls 1
set label "23" at 1.750000, 1.299038 center front
set label "22" at 2.000000, 0.866025 center front
set arrow from 1.750000, 1.299038 to 2.000000, 0.866025 nohead ls 1
set label "22" at 2.000000, 0.866025 center front
set label "23" at 1.750000, 1.299038 center front
set arrow from 2.000000, 0.866025 to 1.750000, 1.299038 nohead ls 1
set label "21" at 1.500000, 0.866025 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 1.500000, 0.866025 to 1.000000, 0.866025 nohead ls 1
set label "22" at 2.000000, 0.866025 center front
set label "18" at 2.500000, 0.866025 center front
set arrow from 2.000000, 0.866025 to 2.500000, 0.866025 nohead ls 1
set label "21" at 1.500000, 0.866025 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 1.500000, 0.866025 to 1.250000, 0.433013 nohead ls 1
set label "23" at 1.750000, 1.299038 center front
set label "15" at 2.000000, 1.732051 center front
set arrow from 1.750000, 1.299038 to 2.000000, 1.732051 nohead ls 1
set label "23" at 1.750000, 1.299038 center front
set label "13" at 1.500000, 1.732051 center front
set arrow from 1.750000, 1.299038 to 1.500000, 1.732051 nohead ls 1
set label "22" at 2.000000, 0.866025 center front
set label "14" at 2.250000, 0.433013 center front
set arrow from 2.000000, 0.866025 to 2.250000, 0.433013 nohead ls 1
set label "21" at 1.500000, 0.866025 center front
set label "20" at 0.750000, 1.299038 center front
set arrow from 1.500000, 0.866025 to 0.750000, 1.299038 nohead ls 2
set label "23" at 1.750000, 1.299038 center front
set label "18" at 2.500000, 0.866025 center front
set arrow from 1.750000, 1.299038 to 2.500000, 0.866025 nohead ls 2
set label "23" at 1.750000, 1.299038 center front
set label "19" at 1.000000, 0.866025 center front
set arrow from 1.750000, 1.299038 to 1.000000, 0.866025 nohead ls 2
set label "22" at 2.000000, 0.866025 center front
set label "20" at 2.750000, 1.299038 center front
set arrow from 2.000000, 0.866025 to 2.750000, 1.299038 nohead ls 2
set label "21" at 1.500000, 0.866025 center front
set label "16" at 1.500000, 0.000000 center front
set arrow from 1.500000, 0.866025 to 1.500000, 0.000000 nohead ls 2
set label "22" at 2.000000, 0.866025 center front
set label "15" at 2.000000, 1.732051 center front
set arrow from 2.000000, 0.866025 to 2.000000, 1.732051 nohead ls 2
set label "22" at 2.000000, 0.866025 center front
set label "17" at 1.250000, 0.433013 center front
set arrow from 2.000000, 0.866025 to 1.250000, 0.433013 nohead ls 2
set label "23" at 1.750000, 1.299038 center front
set label "16" at 2.500000, 1.732051 center front
set arrow from 1.750000, 1.299038 to 2.500000, 1.732051 nohead ls 2
set label "23" at 1.750000, 1.299038 center front
set label "12" at 1.000000, 1.732051 center front
set arrow from 1.750000, 1.299038 to 1.000000, 1.732051 nohead ls 2
set label "21" at 1.500000, 0.866025 center front
set label "14" at 2.250000, 0.433013 center front
set arrow from 1.500000, 0.866025 to 2.250000, 0.433013 nohead ls 2
set label "22" at 2.000000, 0.866025 center front
set label "12" at 2.000000, 0.000000 center front
set arrow from 2.000000, 0.866025 to 2.000000, 0.000000 nohead ls 2
set label "21" at 1.500000, 0.866025 center front
set label "13" at 1.500000, 1.732051 center front
set arrow from 1.500000, 0.866025 to 1.500000, 1.732051 nohead ls 2
plot '-' w d lc 7
0.0 0.0
end
pause -1
