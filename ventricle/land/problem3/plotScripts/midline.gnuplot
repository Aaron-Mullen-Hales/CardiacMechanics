set term pdfcairo dashed enhanced
set datafile commentschars "#"

set size ratio -1
set grid
set xrange [-10:0]
set yrange [-20:5]
set xtics 5
set ytics
set xlabel "x (in mm)"
set ylabel "z (in mm)"
set key off

midline = "midLineDeformed.txt"
reference = "< sed '1d; s/,/ /g' referenceData/land3PlotData.csv"
set style line 1 lc rgb "black" dt 2 lw 1.2

set output "ventricleInflation-midline-deformed.pdf"
plot midline u (1e3*$1):(1e3*$3) w l lw 1.2, \
     reference u 1:2 w l ls 1, \
     midline every 9::0 u (1e3*$1):(1e3*$3) w p pt 7 ps 0.7

set xrange [-2:0]
set yrange [-15:-13]
set size ratio 1
set output "ventricleInflation-midline-apex-deformed.pdf"
plot midline u (1e3*$1):(1e3*$3) w lp pt 7 ps 0.7 lw 1.2

set xrange [-8.75:-8.25]
set yrange [-2:2]
set xtics 0.5
set size ratio 1
set output "ventricleInflation-midline-inflection-deformed.pdf"
plot midline u (1e3*$1):(1e3*$3) w lp pt 7 ps 0.7 lw 1.2
