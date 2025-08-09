# Gnuplot script for xy MCMC results
set terminal pngcairo size 1000,800 enhanced font 'Arial,12'
set output 'figs/xy2d_mh_sweep.png'

set multiplot layout 3,1 title "XY 2D MCMC Beta Sweep\n" font ',14'

# --- 1. Energy by Beta ---
set lmargin 10
set rmargin 5
set tmargin 2
set bmargin 2
set xlabel 'Beta'
set ylabel 'Energy'
set grid
set key outside
plot 'data/xy2d_mh_sweep.dat' using 1:2 with lines title 'Energy'

# --- 2. Magnetization by Beta ---
set xlabel 'Beta'
set ylabel 'Helicity Modulus'
set grid
set key outside
plot 'data/xy2d_mh_sweep.dat' using 1:3 with lines title 'Helicity Modulus'

# --- 3. Autocorrelation Time by Beta ---
set bmargin 5
set xlabel 'Beta'
set ylabel 'Autocorrelation Time'
set grid
set key outside
plot 'data/xy2d_mh_sweep.dat' using 1:4 with lines title 'Helicity Modulus Autocorrelation Time'
