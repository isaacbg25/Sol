set terminal pngcairo size 1000,600 font "Helvetica, 18"
set output "equinocis_esferoide_comp.png"

set xlabel "angle nu (est-oest)"
set ylabel "angle eta (vertical)"

set xrange[-120:120]

set title 'Moviment sol pels diferents equinocis (Terra com a esferoide)'

plot "equinoci1.txt" w l lc 'blue' lw 7 title 'hivern esfera', \
     "equinoci2.txt" w l lc 'red' lw 7 title 'primavera esfera', \
     "equinoci3.txt" w l lc 'yellow' lw 7 title 'estiu esfera', \
     "equinoci4.txt" w l lc 'green' lw 7 title 'tardor esfera', \
     "eqesferoide1.txt" w l lc 'cyan' lw 3 title 'hivern esferoide', \
     "eqesferoide2.txt" w l lc 'pink' lw 3 title 'primavera esferoide', \
     "eqesferoide3.txt" w l lc 'orange' lw 3 title 'estiu esferoide', \
     "eqesferoide4.txt" w l lc 'purple' lw 3 title 'tardor esferoide'
     