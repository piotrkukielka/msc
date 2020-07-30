set palette
set xlabel "time [h]"
set ylabel "distance [km]"
set title "O2 concentration"
set view 50, 300, 1, 1
splot 'results/-1.000000_0.030000' u 1:2:3 w pm3d title ""