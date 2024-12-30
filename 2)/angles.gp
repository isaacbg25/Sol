set terminal wxt
set xlabel "angle nu (est-oest)"
set ylabel "angle eta (vertical)"
set xrange[-120:120]
plot "equinoci1.txt" w l lc 'blue' lw 2 title 'solstici hivern', "equinoci2.txt" w l lc 'red' lw 2 title 'equinoci primavera', "equinoci3.txt" w l lc 'yellow' lw 2 title 'solstici estiu', "equinoci4.txt" w l lc 'green' lw 2 title 'equinoci tardor'
set terminal png
set output "equinocis.png"
replot
set output