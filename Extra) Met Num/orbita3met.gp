set terminal pngcairo size 1000,600 font "Helvetica, 18"

set style line 1 lc rgb "red" #euler
set style line 2 lc rgb "green" #RK2
set style line 3 lc rgb "blue" #RK4
set style line 4 lc rgb "yellow" #sol

set output 'orbitaEul.png'
set ylabel 'y(m)'
set xlabel 'x(m)'
set xrange[-2e11:2e11]
set yrange[-2e11:2e11]
set grid
set title 'Òrbita terrestre (Euler)'
plot 'orbita.txt' using 1:2 with lines lc 'red' lw 2 title 'orbita','sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle

set output 'orbitaRK2.png'
set ylabel 'y(m)'
set xlabel 'x(m)'
set xrange[-2e11:2e11]
set yrange[-2e11:2e11]
set grid
set title 'Òrbita terrestre (RK-2)'
plot 'orbitaRK2.txt' using 1:2 with lines lc 'green' lw 2 title 'orbita','sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle

set output 'orbitaRK4.png'
set ylabel 'y(m)'
set xlabel 'x(m)'
set xrange[-2e11:2e11]
set yrange[-2e11:2e11]
set grid
set title 'Òrbita terrestre (RK-4)'
plot 'orbitaRK4.txt' using 1:2 with lines lc 'blue' lw 2 title 'orbita','sol.txt' using 1:2:3 with labels offset 1,1 notitle,'sol.txt' using 1:2 with points lc rgb "yellow" pt 7 ps 1.5 notitle

set output 'orbita3met.png'
set ylabel 'y(m)'
set xlabel 'x(m)'
set xrange[-2e11:2e11]
set yrange[-2e11:2e11]
set grid
set title 'Òrbita terrestre (comparació 3 mètodes numèrics)'
plot 'orbita.txt' using 1:2 with lines linestyle 1 lw 8 title 'Euler', \
     'orbitaRK2.txt' using 1:2 with lines linestyle 2 lw 5 title 'RK-2', \
     'orbitaRK4.txt' using 1:2 with lines linestyle 3 lw 2 title 'RK-4' , \
     'sol.txt' using 1:2:3 with labels offset 1,1 notitle, 'sol.txt' using 1:2 with points linestyle 4 pt 7 ps 1.5 notitle
