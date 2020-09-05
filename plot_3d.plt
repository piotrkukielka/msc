set palette

# set yrange [3:3.3]
# set xrange[1:1.2]

set xlabel "time [h]"
set ylabel "distance [km]"
set title "DO concentration"
set view 50, 300, 1, 1
splot 'results/-1.000000_0.030000' u 1:2:3 w pm3d title ""
