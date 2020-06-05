set macro  # Enable macro definition

# Some macro definitions

# Colors: hex #RRGGBB
label_color = "#867961"
tic_color = "#383838"
title_color = "#383838"
myblue = "#2196f3"
myred = "#f44336"
mygreen = "#2e7d32"
mypurple = "#9c27b0"
myviolet = "#673ab7"
mybrown = "#795548"
myorange = "#ff9800"

# Width and height of postscript figure in inches
width = 5
height = 8

# Line styles
#set style line 1 lt rgb myblue_color lw 1   # Define linestyle 1
#set style line 2 lt rgb myred_color lw 1    # Define linestyle 2
#set style line 3 lt rgb mygreen_color lw 1  # Define linestyle 3

# Point styles
set style line 1 lc rgb myorange pt 5   # square
set style line 2 lc rgb "purple" pt 7   # circle
set style line 3 lc rgb 'dark-orange' pt 9   # triangle
set style line 4 lc rgb "red" pt 7   # circle
set style line 5 lc rgb 'blue' pt 7   # circle
set style line 6 lc rgb myblue pt 7   # circle
set style line 7 lc rgb myred pt 7   # circle
set style line 8 lc rgb mygreen pt 7   # circle
set style line 9 lc rgb 'black' pt 7   # circle
set style line 10 lc rgb 'gray' pt 7   # circle
set style line 11 lc rgb mybrown pt 7   # circle
set style line 12 lc rgb 'violet' pt 7   # circle
set style line 13 lc rgb 'magenta' pt 7   # circle
set style line 14 lc rgb 'brown' pt 7   # circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set title 'Mass fractions' textcolor rgb title_color
set xlabel 'Log time' textcolor rgb tic_color
set ylabel 'X' textcolor rgb tic_color

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set xrange [-16:-1]
set yrange[-14:1]

set pointsize 0.5    # Size of the plotted points

set key outside right   # place legend

#unset key            # Don't show legend
#set timestamp       # Date/time

# Read data from file and plot to screen

data = "gnuPlot2.out"

plot data using 1:2 ls 1 with lines title 'dt'
replot data using 1:11 ls 2 with lines title '4He'
replot data using 1:12 ls 3 with lines title '12C'
replot data using 1:13 ls 4 with lines title '16O'
replot data using 1:14 ls 5 with lines title '20Ne'
replot data using 1:15 ls 6 with lines title '24Mg'
replot data using 1:16 ls 7 with lines title '28Si'
replot data using 1:17 ls 8 with lines title '32S'
replot data using 1:18 ls 9 with lines title '36Ar'
replot data using 1:19 ls 10 with lines title '40Ca'
replot data using 1:19 ls 11 with lines title '44Ti'
replot data using 1:20 ls 12 with lines title '48Cr'
replot data using 1:21 ls 13 with lines title '52Fe'
replot data using 1:22 ls 14 with lines title '56Ni'

data2 = "gnuPE.inp"

replot data2 using 1:2 ls 1 with lines title 'dt'
replot data2 using 1:11 ls 2 with points title '4He'
replot data2 using 1:12 ls 3 with points title '12C'
replot data2 using 1:13 ls 4 with points title '16O'
replot data2 using 1:14 ls 5 with points title '20Ne'
replot data2 using 1:15 ls 6 with points title '24Mg'
replot data2 using 1:16 ls 7 with points title '28Si'
replot data2 using 1:17 ls 8 with points title '32S'
replot data2 using 1:18 ls 9 with points title '36Ar'
replot data2 using 1:19 ls 10 with points title '40Ca'
replot data2 using 1:19 ls 11 with points title '44Ti'
replot data2 using 1:20 ls 12 with points title '48Cr'
replot data2 using 1:21 ls 13 with points title '52Fe'
replot data2 using 1:22 ls 14 with points title '56Ni'

#plot "scaling_isotopes.dat" using 1:7 ls 2 with points title 'ms/step'
#replot "scaling_isotopes.dat" using 1:7 ls 2 with lines  # fermi gpu

#replot "scaling_isotopes.dat" using 1:6 ls 5 with points  # cpu implicit
#replot "scaling_isotopes.dat" using 1:6 ls 5 with lines  # cpu implicit
 
# Plot to postscript file

set out "massFraction.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 32
replot               # Plot to postscript file

# Plot to PNG file

set out "massFraction.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal pngcairo transparent size 2*width*72, 2*height*72 lw 2 enhanced font 'Arial,28'
replot

quit