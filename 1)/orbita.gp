set terminal png
set output 'orbita.png' 
set ylabel 'y(m)'
set xlabel 'x(m)'
set grid
set title 'Òrbita terrestre'
plot 'orbita.txt' using 1:2 with lines lc 'red' lw 2 title 'orbita','afeli.txt' using 1:2:3 with labels offset 1,1 notitle,'afeli.txt' using 1:2 with points lc rgb "blue" pt 7 ps 1.5 notitle,'periheli.txt' using 1:2:3 with labels offset 1,1 notitle ,'periheli.txt' using 1:2 with points lc rgb "blue" pt 7 ps 1.5 notitle,'sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle
