# Gnuplot script for Ising MCMC results
set terminal pngcairo size 1000,800 enhanced font 'Arial,12'
set output 'figs/ising2d_mh_summary.png'

set multiplot layout 3,1 title "Ising 2D MCMC Results\n" font ',14'

# --- 1. Energy time series ---
set lmargin 10
set rmargin 5
set tmargin 2
set bmargin 2
set xlabel 'Sample Index'
set ylabel 'Energy'
set grid
set key outside
plot 'data/ising2d_mh.dat' using 1:2 with lines title 'Energy'

# --- 2. Magnetization time series ---
set xlabel 'Sample Index'
set ylabel 'Magnetization'
set grid
set key outside
plot 'data/ising2d_mh.dat' using 1:3 with lines title 'Magnetization'

# --- 3. Autocorrelation function ---
set bmargin 5
set xlabel 'Lag / MCS'
set ylabel 'Autocorrelation'
unset y2tics
plot 'data/ising2d_mh_autocorr.dat' using 1:2 with lines title 'Magnetization Autocorr.'

