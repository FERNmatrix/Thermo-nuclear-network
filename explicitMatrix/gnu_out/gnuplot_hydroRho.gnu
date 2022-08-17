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
width = 10.5
height = 6.5

set pointsize 1.5    # Size of the plotted points

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
set style line 6 lc rgb 'red'  pt 6   # circle
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

set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'log density (g/cm^3)' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

set timestamp       # Date/time

ds="tidal SN profile"
ds = ds.": Density"
set title ds textcolor rgb title_color


# -------- Axis ranges and ticmarks -----------

xlow = -3
xup = 6
xtics = 1    # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = 0
yup = 1e4
ytics = 1000    # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

# set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


#file1 = "gnufile2.data"
file1 = "plot4.data"

#plot file1 using 1:10 with lines ls 1 lw 1.0 dashtype 1 title "T interp (t)"
#plot file1 using 1:(10**$11) with lines ls 1 lw 1.0 dashtype 1 title "Rho interp (t)"
#plot file1 using 1:11 with lines ls 1 lw 1.0 dashtype 1 title "log Rho interp (t)"
plot file1 using 1:(10**$3) with lines ls 1 lw 1.0 dashtype 1 title "Rho interp (t)"

file2 = "hydroProfile.out"  # The input hydro profile

#replot file2 using (log10($1)):(log10($2)) with points ls 11 lw 0.5 dashtype 1 title "T input (t)"
#replot file2 using (log10($1)):($2/1e9) with points ls 11 lw 0.5 dashtype 1 title "T input (t)"
#replot file2 using 1:2 with lines ls 11 lw 0.5 dashtype 2 title "T input (t)"
#replot file2 using 1:(10**$3) with points ls 6 lw 0.5 dashtype 2 title "Rho input (t)"
#replot file2 using 1:3 with points ls 6 lw 1.0 dashtype 2 title "log Rho input (t)"
replot file2 using 1:(10**$3) with points ls 6 lw 1.0 dashtype 2 title "Rho input (t)"


# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,20"
set key top right font "Arial,20"
set xlabel 'log time' textcolor rgb tic_color font "Arial,22"
set ylabel 'log density (g/cm^3)' textcolor rgb tic_color font "Arial,22"

# Plot to postscript file

set out "gnuplot_hydroRho.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 #"Arial" 24
replot               # Plot to postscript file

# Plot to PNG file

#set out "gnuplot_hydroRho.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit
