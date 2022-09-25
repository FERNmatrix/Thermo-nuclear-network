set macro  # Enable macro definition

# Some macro definitions

# Colors: hex #RRGGBB
label_color = "#000000"
tic_color = "#000000"
title_color = "#000000"
myblue = "#2196f3"
myred = "#f44336"
mygreen = "#2e7d32"
mypurple = "#9c27b0"
myviolet = "#673ab7"
mybrown = "#795548"
myorange = "#ff9800"

# Width and height of postscript figure in inches
width = 8
height = 8

# x-axis resolution
set samples 1000

# Line styles.  
# For lines: plot x with lines ls 1
# For points: plot x with points ls 1

set style line 1 lc rgb 'black' pt 5   # fill square
set style line 2 lc rgb myred pt 7   # circle
set style line 3 lc rgb 'blue' pt 9   # triangle
set style line 4 lc rgb mygreen pt 7   # circle
set style line 5 lc rgb mypurple pt 7   # circle
set style line 6 lc rgb myviolet  pt 7   # circle
set style line 7 lc rgb mybrown  pt 7   # circle
set style line 8 lc rgb myorange  pt 7   # circle
set style line 9 lc rgb myblue pt 7   # circle
set style line 10 lc rgb 'green' pt 6   # open circle
set style line 11 lc rgb 'red' pt 4   # open circle
set style line 12 lc rgb 'gray' pt 4   # open circle
set style line 13 lc rgb 'gold' pt 4   # open circle
set style line 14 lc rgb 'orange' pt 4   # open circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log t (s)' textcolor rgb tic_color #font "Arial,32"
set ylabel 'Log dt (s)' textcolor rgb tic_color #font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key outside    # Place legend outside top
#unset key            # Don't show legend

set timestamp       # Date/time

ds="C++ Asy and Asy+PE alpha with PF"
ds = ds.": T9=5 rho=1e8"
set title ds textcolor rgb title_color



# -------- Axis ranges and ticmarks -----------

xlow = -16
xup = 2
xtics = 1    # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -17 
yup = 2
ytics = 2      # Space between major y ticmarks
minytics = 5   # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------

file1 = "plot2.data"     # C++
#file2 = "dataRef/gnufile_alpha_T9_5_1e8_asy+pe_java.data"  # Java reference
file3 = "dataRef/gnufile_alpha_T9_5_1e7_asy.data"   # C++ asy

# Plot C++ Asy+PE data

plot file1 using 1:5 with lines ls 2 lw 1.0 dashtype 2 title "log 2/Rmax"
#replot file1 using 1:($7+0.301) with lines ls 4 lw 1.0 dashtype 2 title "log 2*dt-FE"
replot file1 using 1:2 with lines ls 3 lw 2.0 dashtype 1 title "log dt (C++ asy+pe)"

# Plot Java Asy+PE data

#replot file2 using 1:2 with lines ls 8 lw 2.0 dashtype 1 title "log dt (java asy+pe)"

# Plot C++ asy

replot file3 using 1:2 with lines ls 10 lw 2.0 dashtype 1 title "log dt (c++ asy)"

# Plot dt contours

replot file1 using 1:( log10((10**$1)*0.1) ) with lines ls 1 lw 1.0 dashtype 2 title "dt=0.1 t"
replot file1 using 1:( log10((10**$1)*0.01) ) with lines ls 1 lw 1.0 dashtype 0 title "dt=0.01 t"
replot file1 using 1:( log10((10**$1)*0.001) ) with lines ls 1 lw 1.0 dashtype 7 title "dt=0.001 t"
replot file1 using 1:( log10((10**$1)*0.0001) ) with lines ls 1 lw 1.0 dashtype 9 title "dt=0.0001 t"
#replot file1 using 1:( log10((10**$1)*0.00001) ) with lines ls 1 lw 1.0 dashtype 8 title "dt=0.00001 t"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top left inside font "Arial,16"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,22"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,22"

# Plot to postscript file

set out "gnuplot_stiffstuff.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Symbol" 22
replot               # Plot to postscript file

# Plot to PNG file

#set out "gnuplot_stiffstuff.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit
