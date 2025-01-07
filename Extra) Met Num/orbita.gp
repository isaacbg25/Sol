set terminal png
set output 'orbitaE.png' 
set ylabel 'y(m)'
set xlabel 'x(m)'
set grid
set xrange[-2e11:2e11]
set yrange[-2e11:2e11]
set title 'Òrbita (Euler)'
plot 'orbita.txt' using 1:2 with lines lc 'red' lw 2 title 'orbita','sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle
