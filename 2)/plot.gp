set terminal pngcairo size 1000, 600 font "Helvetica, 18"

set output "equinocis.png"
set xlabel "Angle horitzontal (est-oest) (°)"
set ylabel "Angle vertical (°)"
set xrange[-120:120]
plot "equinoci1.txt" w l lc 'blue' lw 2 title 'Solstici hivern', "equinoci2.txt" w l lc 'red' lw 2 title 'Equinoci primavera', "equinoci3.txt" w l lc 'yellow' lw 2 title 'Solstici estiu', "equinoci4.txt" w l lc 'green' lw 2 title 'Equinoci tardor'