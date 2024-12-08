set terminal png
set output 'orbita.png' 

plot 'orbita.txt' using 1:2 with lines title 'orbita',\
     'punts.txt' using 1:2 with points pointtype 7 pointsize 1.5 linecolor rgb "black", \
     'punts.txt' using 1:2:3 with labels textcolor variable font ",12" offset 0.5,0.5 notitle


