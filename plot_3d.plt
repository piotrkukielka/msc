set terminal png
set output '2_zoom.png'
set encoding utf8

set lmargin 0.5
set rmargin 1
set tmargin 0
set bmargin 0

set cblabel "stężenie"

set ztics -1,0.2,2

set palette

set xlabel "czas" offset 2, 0, 0
set ylabel "odległość" offset 0, -0.5, 0
set zlabel "stężenie" rotate parallel offset 2, 0, 0

# set view 50, 300, 1, 1
# set view 50, 180, 1, 1 # behind view
#splot 'results/-1.000000_0.030000' u 1:2:3 w pm3d title ""

# zoom
set yrange [3:3.15]
set xrange[1:1.2]
set view 80, 275, 1, 1
set xtics -1,0.1,2
set ytics offset 0,-1
set ylabel offset 0, -2, 0
splot 'real_results/old_preparameters/peclet/pe2' u 1:2:3 w pm3d title ""