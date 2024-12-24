set terminal wxt
set xlabel "angle est-oest"
set ylabel "angle vertical"
plot "angles1.txt" w l lc 'blue', "angles2.txt" w l lc 'red', "angles3.txt" w l lc 'yellow'
set terminal png
set output "anglesvaris.png"
replot
set output