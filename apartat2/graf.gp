set title "3D Points from data.txt"
set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Z-axis"
set grid
splot "data.txt" using 1:2:3 with lines  title "prova"
set terminal pngcairo
set output "prova2.png"
replot
set output
