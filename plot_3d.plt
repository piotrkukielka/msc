set palette

set xlabel "time [h]"
set ylabel "distance [km]"
set title "DO concentration"
set view 50, 300, 1, 1
splot 'results/-1.000000_0.030000' u 1:2:3 w pm3d title ""

# zoom
#set yrange [3:3.15]
#set xrange[1:1.2]
#set view 80, 275, 1, 1
#splot 'real_results/peclet/pe2.4' u 1:2:3 w pm3d title ""