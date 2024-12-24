set terminal wxt
set xlabel "angle est-oest"
set ylabel "angle vertical"
plot "equinoci1.txt" w l lc 'blue' lw 2 title 'equinoci primavera', "equinoci2.txt" w l lc 'red' lw 2 title 'solstici estiu', "equinoci3.txt" w l lc 'yellow' lw 2 title 'equinoci tardor', "equinoci4.txt" w l lc 'green' lw 2 title 'solstici hivern'
set terminal png
set output "equinocis.png"
replot
set output