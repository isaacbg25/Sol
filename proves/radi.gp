set terminal png
set output 'radi.png' 

plot 'radi.txt' using 1:2 with lines title 'radi(t)'