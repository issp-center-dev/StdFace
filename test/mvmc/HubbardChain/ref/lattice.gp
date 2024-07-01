#set terminal pdf color enhanced \
#dashed dl 1.0 size 20.0cm, 20.0cm 
#set output "lattice.pdf"
set xrange [-2.000000: 8.000000]
set yrange [-2.000000: 8.000000]
set size square
unset key
unset tics
unset border
set style line 1 lc 1 lt 1
set style line 2 lc 5 lt 1
set style line 3 lc 0 lt 1
set arrow from 0.000000, 0.000000 to 1.000000, 0.000000 nohead front ls 3
set arrow from 1.000000, 0.000000 to 1.000000, 6.000000 nohead front ls 3
set arrow from 1.000000, 6.000000 to 0.000000, 6.000000 nohead front ls 3
set arrow from 0.000000, 6.000000 to 0.000000, 0.000000 nohead front ls 3
set label "0" at 0.000000, 0.000000 center front
set label "5" at 0.000000, -1.000000 center front
set arrow from 0.000000, 0.000000 to 0.000000, -1.000000 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "1" at 0.000000, 1.000000 center front
set arrow from 0.000000, 0.000000 to 0.000000, 1.000000 nohead ls 1
set label "0" at 0.000000, 0.000000 center front
set label "4" at 0.000000, -2.000000 center front
set arrow from 0.000000, 0.000000 to 0.000000, -2.000000 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "2" at 0.000000, 2.000000 center front
set arrow from 0.000000, 0.000000 to 0.000000, 2.000000 nohead ls 2
set label "0" at 0.000000, 0.000000 center front
set label "3" at 0.000000, -3.000000 center front
set label "0" at 0.000000, 0.000000 center front
set label "3" at 0.000000, 3.000000 center front
set label "1" at 0.000000, 1.000000 center front
set label "0" at 0.000000, 0.000000 center front
set arrow from 0.000000, 1.000000 to 0.000000, 0.000000 nohead ls 1
set label "1" at 0.000000, 1.000000 center front
set label "2" at 0.000000, 2.000000 center front
set arrow from 0.000000, 1.000000 to 0.000000, 2.000000 nohead ls 1
set label "1" at 0.000000, 1.000000 center front
set label "5" at 0.000000, -1.000000 center front
set arrow from 0.000000, 1.000000 to 0.000000, -1.000000 nohead ls 2
set label "1" at 0.000000, 1.000000 center front
set label "3" at 0.000000, 3.000000 center front
set arrow from 0.000000, 1.000000 to 0.000000, 3.000000 nohead ls 2
set label "1" at 0.000000, 1.000000 center front
set label "4" at 0.000000, -2.000000 center front
set label "1" at 0.000000, 1.000000 center front
set label "4" at 0.000000, 4.000000 center front
set label "2" at 0.000000, 2.000000 center front
set label "1" at 0.000000, 1.000000 center front
set arrow from 0.000000, 2.000000 to 0.000000, 1.000000 nohead ls 1
set label "2" at 0.000000, 2.000000 center front
set label "3" at 0.000000, 3.000000 center front
set arrow from 0.000000, 2.000000 to 0.000000, 3.000000 nohead ls 1
set label "2" at 0.000000, 2.000000 center front
set label "0" at 0.000000, 0.000000 center front
set arrow from 0.000000, 2.000000 to 0.000000, 0.000000 nohead ls 2
set label "2" at 0.000000, 2.000000 center front
set label "4" at 0.000000, 4.000000 center front
set arrow from 0.000000, 2.000000 to 0.000000, 4.000000 nohead ls 2
set label "2" at 0.000000, 2.000000 center front
set label "5" at 0.000000, -1.000000 center front
set label "2" at 0.000000, 2.000000 center front
set label "5" at 0.000000, 5.000000 center front
set label "3" at 0.000000, 3.000000 center front
set label "2" at 0.000000, 2.000000 center front
set arrow from 0.000000, 3.000000 to 0.000000, 2.000000 nohead ls 1
set label "3" at 0.000000, 3.000000 center front
set label "4" at 0.000000, 4.000000 center front
set arrow from 0.000000, 3.000000 to 0.000000, 4.000000 nohead ls 1
set label "3" at 0.000000, 3.000000 center front
set label "1" at 0.000000, 1.000000 center front
set arrow from 0.000000, 3.000000 to 0.000000, 1.000000 nohead ls 2
set label "3" at 0.000000, 3.000000 center front
set label "5" at 0.000000, 5.000000 center front
set arrow from 0.000000, 3.000000 to 0.000000, 5.000000 nohead ls 2
set label "3" at 0.000000, 3.000000 center front
set label "0" at 0.000000, 0.000000 center front
set label "3" at 0.000000, 3.000000 center front
set label "0" at 0.000000, 6.000000 center front
set label "4" at 0.000000, 4.000000 center front
set label "3" at 0.000000, 3.000000 center front
set arrow from 0.000000, 4.000000 to 0.000000, 3.000000 nohead ls 1
set label "4" at 0.000000, 4.000000 center front
set label "5" at 0.000000, 5.000000 center front
set arrow from 0.000000, 4.000000 to 0.000000, 5.000000 nohead ls 1
set label "4" at 0.000000, 4.000000 center front
set label "2" at 0.000000, 2.000000 center front
set arrow from 0.000000, 4.000000 to 0.000000, 2.000000 nohead ls 2
set label "4" at 0.000000, 4.000000 center front
set label "0" at 0.000000, 6.000000 center front
set arrow from 0.000000, 4.000000 to 0.000000, 6.000000 nohead ls 2
set label "4" at 0.000000, 4.000000 center front
set label "1" at 0.000000, 1.000000 center front
set label "4" at 0.000000, 4.000000 center front
set label "1" at 0.000000, 7.000000 center front
set label "5" at 0.000000, 5.000000 center front
set label "4" at 0.000000, 4.000000 center front
set arrow from 0.000000, 5.000000 to 0.000000, 4.000000 nohead ls 1
set label "5" at 0.000000, 5.000000 center front
set label "0" at 0.000000, 6.000000 center front
set arrow from 0.000000, 5.000000 to 0.000000, 6.000000 nohead ls 1
set label "5" at 0.000000, 5.000000 center front
set label "3" at 0.000000, 3.000000 center front
set arrow from 0.000000, 5.000000 to 0.000000, 3.000000 nohead ls 2
set label "5" at 0.000000, 5.000000 center front
set label "1" at 0.000000, 7.000000 center front
set arrow from 0.000000, 5.000000 to 0.000000, 7.000000 nohead ls 2
set label "5" at 0.000000, 5.000000 center front
set label "2" at 0.000000, 2.000000 center front
set label "5" at 0.000000, 5.000000 center front
set label "2" at 0.000000, 8.000000 center front
plot '-' w d lc 7
0.0 0.0
end
pause -1
