set terminal png
set output 'orbita.png' 

plot 'orbita.txt' using 1:2 with lines title 'orbita',\
     '-' using 1:2 with points pointtype 7 pointsize 1.5 linecolor rgb "black", \
     '' using 1:2:3:4 with labels textcolor variable font ",12" offset 0.5,0.5 notitle

0 0 "sol" "yellow"
-1.47e11 0 "afeli" "blue"
1.52e11 0 "periheli" "blue"

