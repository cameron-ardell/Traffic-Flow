set xlabel 'x'
set ylabel 't'
set title 'Rho vs x'
set xrange[0:100]
set yrange[0:1]
splot "Rho over x and t" u 1:2:3 w l,
