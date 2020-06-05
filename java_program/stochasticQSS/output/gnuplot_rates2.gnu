set macro  # Enable macro definition

# Some macro definitions

# Custom colors: hex #RRGGBB

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
set style line 14 lc rgb 'dark-grey' pt 4   # open circle
set style line 15 lc rgb 'dark-cyan' pt 4   # open circle
set style line 16 lc rgb 'purple' pt 4   # open circle
set style line 17 lc rgb 'orchid' pt 4   # open circle
set style line 18 lc rgb 'brown' pt 4   # open circle
set style line 19 lc rgb 'dark-green' pt 4   # open circle
set style line 20 lc rgb 'magenta' pt 4   # open circle

# Width and height of postscript figure in inches

width = 6.0 
height = 8.5 

#width = 8.5
#height = 5.0

# x-axis plot resolution
set samples 2000

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log T (K)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log Rate' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
set logscale x
set logscale y

set pointsize 1.5    # Size of any plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

set timestamp       # Date/time

# Input data file
file1 = "gnuplotRates2.data"

# Main title
ds="Rates: "

# Axis ranges and ticmarks
set yrange[1e-20 : 1e10]
set ytics 1e-20, 10, 1e10
#set mytics 10 # set minor y tics
set xrange [1e8 : 1e10]
#set xtics  1e8, 10, 1e10
#set mxtics 5 # set minor x tics

set grid   # set grid at major ticmarks

#####################################################################################
# Read in data from file and assign to gnuplot variables.  Different
# technique for string variables and floating point variables.  Both
# methods strange because gnuplot only knows how to plot file input data, so
# it has to be fooled into using the plot command to store data in a file
# row and column in gnuplot number or string variable.
#####################################################################################

#####################################################################################
# Reading in a string variable.  General form is
#####################################################################################
#
#   set terminal unknown  # No plot will be generated
#   plot filename every ::textRow-1::textRow-1 using (textVar=stringcolumn(textColumn),0)
#
# where filename is the file, and the string desired is in row textRow and column
# textColumn of the data file, which will be assigned to gnuPlot variable textVar. 
# After setting terminal unknown, must then reset to a plotting terminal to get plots.
# CAUTION: In gnuplot row numbers start with 0 but column numbers start with 1 (because
# column 0 is the line number?).

######################################################################################################
# Reading in numerical variables from file. 
######################################################################################################
# See https://stackoverflow.com/questions/11211339/gnuplot-store-one-number-from-data-file-into-variable
# Gnuplot can access the data only through the "plot" command:
# Syntax:
#
# set table
# plot file u 0:($0==RowIndex?(VariableName1=$ColumnIndex):$ColumnIndex),\
#        '' u 0:($0==RowIndex?(VariableName2=$ColumnIndex):$ColumnIndex),\
#                .
#                .
#                .
# unset table
#
# where RowIndex starts with 0, ColumnIndex starts with 1,
# 'u' is an abbreviation for the 'using' modifier, and file
# is the filename containing the data.
#
# EXAMPLE:
#
# Loop over lines in data file and extract the pn[i] in its rows
# numberCases=3
#set table
#do for [i=1:numberCases]{
#j=i-1
#plot file1 u 0:($0==j?(p1[i]=$2):$2),\
#        '' u 0:($0==j?(p2[i]=$3):$3),\
#        '' u 0:($0==j?(p3[i]=$4):$4),\
#        '' u 0:($0==j?(p4[i]=$5):$5),\
#        '' u 0:($0==j?(p5[i]=$6):$6),\
#        '' u 0:($0==j?(p6[i]=$7):$7),\
#        '' u 0:($0==j?(p7[i]=$8):$8),\
#        '' u 0:($0==j?(reaction[i]=$1):$1)
#}
#unset table
######################################################################################################


#####-- Read in data from file --#####

# Read in some numerical variables from file for the number
# of reactions numReactions and number of plot points for each curve
# numPoints. gnuplot seems to think these are floats, but they work
# fine as loop indices.

set table
plot file1 u 0:($0==0?(numReactions=$2):$2),\
     '' u 0:($0==0?(numPoints=$3):$3)
unset table

# Read in the reaclib classes for each reaction
array reaclibClass[60]

#plot file u 0:($0==RowIndex?(VariableName1=$ColumnIndex):$ColumnIndex),\

#plot file1 u 0:($0==1?(reaclibClass[1]=$1):$1,\
#        '' u 0:($0==1?(reaclibClass[2]=$2):$2)

# Not sure how to get this to work in a loop, so unroll
set table
plot file1 u 0:($0==1?(reaclibClass[1]=$1):$1),\
     '' u 0:($0==1?(reaclibClass[2]=$2):$2)

unset table


# Create an array to hold the reaction strings for each reaction. Although the
# array reaction[] will be read in with a stringColumn method in gnuplot, unless
# the array is initialized as a string, gnuplot will have trouble realizing that it
# is a string in the plotting loop below when used for the title. Hence the loop to
# initialize to an empty string.

array reaction[60]
do for [i=1:20]{
    reaction[i]=""
}

set terminal unknown   # No plot will be generated

# Read in title
plot file1 every ::0::0 using (ratetitle=stringcolumn(1),0)

# Read in reaction labels
do for [i=1:numReactions]{
    plot file1 every ::2::2 using (reaction[i]=stringcolumn(i+1),0)
}

set terminal x11   # Reset to plotting terminal for later plots

unitString=" units: s^{-1}(1-body); cm^3mol^{-1}s^{-1}(2-body); cm^6mol^{-2}s^{-1}(3-body)"
ds = ratetitle.unitString
set title ds textcolor rgb title_color

####--Plot the curves--####

# Use a plotting loop to plot all curves.  A * argument for
# the loop means to process until no more lines in data file.
# Note the initialization of reaction[] above to an empty string.
# Otherwise gnuplot throws a not a string error in trying to use
# reaction[i-1] as a legend for each curve, though it work fine
# for an explicit plot command such as
#    plot file1 using 1:2 with lines ls 1 title reaction[1]
# without such an initialization.

plot for [i=2:numReactions+1] file1 using 1:i with lines ls i title reaction[i-1]

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,21"
set key top right font "Arial,21"


####--Plot to postscript file--####

set out "gnuplot_rates2.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 26
replot               # Plot to postscript file


####--Plot to PNG file--####

set out "gnuplot_rates2.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2
replot

quit
