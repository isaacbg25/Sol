set terminal pngcairo size 1000,600 font "Helvetica, 18"
set output "equinocis_esferoide.png"

set xlabel "angle nu (est-oest)"
set ylabel "angle eta (vertical)"

set xrange[-120:120]

set title 'Moviment sol pels diferents equinocis (Terra com a esferoide)'

plot "eqesferoide1.txt" w l lc 'blue' lw 2 title 'solstici hivern', \
     "eqesferoide2.txt" w l lc 'red' lw 2 title 'equinoci primavera', \
     "eqesferoide3.txt" w l lc 'yellow' lw 2 title 'solstici estiu', \
     "eqesferoide4.txt" w l lc 'green' lw 2 title 'equinoci tardor'

