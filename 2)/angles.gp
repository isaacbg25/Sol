set terminal wxt
set xlabel "angle est-oest"
set ylabel "angle vertical"
plot "angles1.txt" w p lc 'blue', "angles2.txt" w p lc 'red', "angles3.txt" w p lc 'yellow'
set terminal png
set output "anglesvaris.png"
replot
set output