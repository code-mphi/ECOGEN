# Plot experimental and numerical saturation curve p-v
# Numerical saturation curve can be obtained with the 
# script psat-num.py
# Experimental data can be found in NIST Fluid database

# Be aware that numerical pressure is converted to bar

reset

set key font ",12"
set title font ",12"
set xlabel font ",12"
set tics font ",10"

set style data points
set grid
set key

# Zoom on liquid/gas 
set xrange [-2:102]
set yrange [0:0.12]

# Zoom on gas 
# set xrange [3.6:150]
# set yrange [0:0.12]

# Zoom on liquid
# set xrange [0.0009:0.0011]
# set yrange [0.006:0.1]

set xlabel "v (m3/kg)"
set title "p (bar)"

plot "psat-num.out" u 3:($2*0.00001) lc "blue" lw 3 w p title 'vL - num.', \
  "psat-num.out" u 4:($2*0.00001) lc "forest-green" lw 3 w p title 'vG - num.', \
  "psat-exp.out" u 4:2 title 'vL - exp.' lc "red" lw 1.5 w l, \
  "psat-exp.out" u 16:2 title 'vG - exp.' lc "orange" lw 1.5 w l

pause(-1)
