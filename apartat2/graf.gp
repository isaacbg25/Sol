set title "3D Points from data.txt"
set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Z-axis"
set grid
splot "data.txt" ::1::25 using 1:2:3 with points lc 'blue' title "1",\
      "data.txt" ::27::51 using 1:2:3 with points lc 'red' title "2",\
      'data.txt' ::53::77 using 1:2:3 with points lc 'green' title "3",\
      "data.txt" ::79::103 using 1:2:3 with points lc 'yellow' title "4"
      
#splot "data2.txt" using 1:2:3 with lines lc 'red', "data1.txt" using 1:2:3 with lines lc 'blue', "data3.txt" using 1:2:3 with lines lc 'green'

      
set terminal pngcairo
set output "equador.png"
replot
set output
