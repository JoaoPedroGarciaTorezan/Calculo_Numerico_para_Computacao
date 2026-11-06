set grid
set xlabel 'x'
set ylabel 'y'
plot 'grafico_modo1.dat' u 1:2 w l lw 2 t 'f(x)', '' u 1:3 w l lw 2 dt 2 t 'h(x)'
pause -1
