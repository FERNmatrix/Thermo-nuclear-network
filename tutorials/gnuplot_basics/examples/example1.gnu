set macro  # Enable macro definition

# Some macro definitions

# Custom colors: hex #RRGGBB in quotes

label_color = "#867961"
tic_color = "#383838"
title_color = "#383838"
myblue_color = "#5ea2c6
myred_color = "#bb6255"
mygreen_color = "#668874"

# Width and height of postscript figure in inches
width = 8
height = 5

# Line styles
set style line 1 lt rgb myblue_color lw 1   # Define linestyle 1
set style line 2 lt rgb myred_color lw 1    # Define linestyle 2
set style line 3 lt rgb mygreen_color lw 1  # Define linestyle 3
set style line 4 lt rgb "black" lw 1        # Define linestyle 4
set style line 5 lt rgb "purple" lw 1       # Define linestyle 5
set style line 6 lt rgb "red" lw 1          # Define linestyle 6
set style line 7 lt rgb "royalblue" lw 1    # Define linestyle 7
set style line 8 lt rgb "goldenrod" lw 1    # Define linestyle 8
set style line 9 lt rgb "green" lw 1        # Define linestyle 9
set style line 10 lt rgb "orchid" lw 1       # Define linestyle 10
set style line 11 lt rgb "brown" lw 1       # Define linestyle 11

# Point styles
set style line 12 lc rgb "#5ea2c6" pt 5   # square
set style line 13 lc rgb "#5ea2c6" pt 7   # circle
set style line 14 lc rgb 'dark-orange' pt 9   # triangle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees
# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color
set bmargin 4  # Bottom margin
# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set title 'Landau Dirac Dispersion with B' textcolor rgb title_color
set xlabel 'B' textcolor rgb tic_color
set ylabel 'E' textcolor rgb tic_color

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set xrange [0:1]
set yrange[-10:10]
set pointsize 1.0    # Size of plotted points

#unset key           # Don't show legend
set key out          # Put legend outside plot
set key vert         # vert for vertical; horiz for horizontal
set key top right    # Move legend to upper right

#set timestamp       # Date/time

Omega = 1.0

# Define functions to plot

scale = 1.5
plusE(B,n) = scale*sqrt(n*B)
minusE(B,n) = -scale*sqrt(n*B)
# Make plots

plot[B=0:10] plusE(B,0) ls 1, plusE(B,1) ls 1, plusE(B,2) ls 1, plusE(B,3) ls 1 
replot minusE(B,1) ls 1, minusE(B,2) ls 1, minusE(B,3) ls 1

# Plot to postscript file

set out "example1.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 32
replot               # Plot to postscript file

# Plot to PNG file

set out "example1.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2 enhanced font 'Arial,28'
replot

quit
 

