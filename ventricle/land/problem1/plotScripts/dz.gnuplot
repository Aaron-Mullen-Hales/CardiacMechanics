set term pdfcairo dashed enhanced
set datafile separator " "
set output "land_problem1_dz_convergence.pdf"

set grid
set xrange [1.2:0.01]
set yrange [*:*]
set logscale x
set xtics (1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0138889)

set xlabel "Cell size (in mm)"
set ylabel "Deformed tip z-coordinate (in mm)"
set key top right

# The function object reports displacement.  The monitored point starts at
# z = 0.001 m, so add that initial coordinate before converting to millimetres.
initialTipZ = 0.001
referenceZmm = 4.18

plot \
    "runs/beam.summary.txt" u 2:(1e3*($4 + initialTipZ)) w lp pt 7 lc rgb "#d7191c" t "RhieChow", \
    referenceZmm w l dt 2 lw 2 lc rgb "black" t "Reference = 4.18 mm"
