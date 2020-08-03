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
width = 8.5
height = 4.5

# x-axis resolution.  Larger gives better resolution
# but larger files.

#set samples 1000

# Line styles  

set style line 1 lc rgb myblue pt 5    # fill square
set style line 2 lc rgb myred pt 7     # circle
set style line 3 lc rgb 'purple' pt 9  # triangle

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Theta (deg)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Value' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 0.7    # Size of any plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

set timestamp         # Date/time label

ds="Cosine, sine, and tangent"

file1 = "gnufile.data"    # Input data file

set xrange [0 : 360]
set yrange[-3 : 3]

set title ds textcolor rgb title_color

# Plot data read in from file, with independent variable
# in column 1 and dependent variables in columns 2, 3, 4.

plot file1 using 1:2 with points ls 2 title "Cos"   # Points
replot file1 using 1:3 with lines ls 3 title "Sin"  # Line
replot file1 using 1:4 with lines ls 1 title "Tan"  # Line

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"

# Plot to postscript file

set out "gnuplot_example.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 26
replot               # Plot to postscript file

# Plot to PNG file

set out "gnuplot_example.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2
replot

quit
