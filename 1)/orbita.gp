set terminal pngcairo size 900, 600

set output 'orbita.png' 
set size ratio -1

set ylabel 'y(m)'
set xlabel 'x(m)'
plot 'orbita.txt' using 1:2 with lines lc 'red' lw 2 notitle,'afeli.txt' using 1:2:3 with labels offset 1,1 notitle,'afeli.txt' using 1:2 with points lc rgb "blue" pt 7 ps 1.5 notitle,'periheli.txt' using 1:2:3 with labels offset 1,1 notitle ,'periheli.txt' using 1:2 with points lc rgb "blue" pt 7 ps 1.5 notitle,'sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle
