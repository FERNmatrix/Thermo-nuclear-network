
#**********************************************************PLOT DT vs T****************************************************************#
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
width = 8.5
height = 4.5

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

set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": dt vs t"
set title ds textcolor rgb title_color

file1 = "gnufileBENCH.data"
file2 = "gnufile.data"

# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -18
yup =4
ytics = 1      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------

#---------------------PLOT DT vs T -------------------------#
plot file2 using 1:2 with lines ls 2 title "dt"
#-----------------------------------------------------------#

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_dt.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file


#**********************************************************************************************************************************************************#


#***************************************************************** PLOT FRAC EQUIL ************************************************************************#

set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Frac Asy, Frac RG equil' textcolor rgb tic_color font "Arial,32"

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": Fraction asymptotic and equilibrated vs t"

set title ds textcolor rgb title_color

file3 = "gnufile.data"

# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = 0.0 #-0.01
yup = 1.0 #1.01
ytics = 0.1      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

set grid   # set x-y grid at major ticmarks

# --------------------- PLOT FRAC EQUIL --------------------#
plot file3 using 1:5 with lines ls 2 title "Frac Asy"
replot file3 using 1:6 with lines ls 3 title "Frac RG equil"
#-----------------------------------------------------------#

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Frac Asy, Frac RG equil' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_frac.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

#*****************************************************************************************************************************************************************#

#*************************************************************** PLOT X COMPARISONS ******************************************************************************#

set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": mass fraction"
set title ds textcolor rgb title_color

file4 = "gnufileBENCH.data"
file5 = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -14
yup = 0
ytics = 1.0      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Edit the following plot commands to correspond to data
# read in from data file

plot file4 using 1:8 with lines ls 2 title "4He"
replot file5 using 1:9 with lines ls 1 title "TS-HE"

replot file4 using 1:9 with lines ls 2 title "12C"
replot file5 using 1:11 with lines ls 3 title "TS-C"

replot file4 using 1:10 with lines ls 2 title "16O"
replot file5 using 1:13 with lines ls 4 title "TS-O"

replot file4 using 1:11 with lines ls 2 title "20Ne"
replot file5 using 1:15 with lines ls 5 title "TS-Ne"

replot file4 using 1:12 with lines ls 2 title "24Mg"
replot file5 using 1:17 with lines ls 6 title "TS-Mg"

replot file4 using 1:13 with lines ls 2 title "28Si"
replot file5 using 1:19 with lines ls 7 title "TS-Si"

replot file4 using 1:14 with lines ls 2 title "32S"
replot file5 using 1:21 with lines ls 8 title "TS-S"

replot file4 using 1:15 with lines ls 2 title "36Ar"
replot file5 using 1:23 with lines ls 9 title "TS-Ar"

replot file4 using 1:16 with lines ls 2 title "40Ca"
replot file5 using 1:25 with lines ls 10 title "TS-Ca"

replot file4 using 1:17 with lines ls 2 title "44Ti"
replot file5 using 1:27 with lines ls 18 title "TS-Ti"

replot file4 using 1:18 with lines ls 2 title "48Cr"
replot file5 using 1:29 with lines ls 5 title "TS-Cr"

replot file4 using 1:19 with lines ls 2 title "52Fe"
replot file5 using 1:31 with lines ls 13 title "TS-Fe"

replot file4 using 1:20 with lines ls 2 title "56Ni"
replot file5 using 1:33 with lines ls 14 title "TS-Ni"

replot file4 using 1:21 with lines ls 2 title "60Zn"
replot file5 using 1:35 with lines ls 1 title "TS-Zn"

replot file4 using 1:22 with lines ls 2 title "64Ge"
replot file5 using 1:37 with lines ls 3 title "TS-Ge"

replot file4 using 1:23 with lines ls 2 title "68Se"
replot file5 using 1:39 with lines ls 4 title "TS-Se"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_X.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

#*****************************************************************************************************************************************************************#

#***************************************************************** SUM X ****************************************************************************************#
# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Sum X' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": Sum X  vs t"
set title ds textcolor rgb title_color

file1 = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = 0.98
yup = 1.02
ytics = 0.01      # Space between major y ticmarks
minytics = 4  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------



# Edit the following plot commands to correspond to data
# read in from data file

plot file1 using 1:7 with lines ls 2 title "sum X"


# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Sum X' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_sumX.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

#************************************************************************************************************************************************************#

#***************************************************************** DT vs RATES ******************************************************************************#
# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": dt vs t"
set title ds textcolor rgb title_color

file1 = "gnufile2.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 4
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -18
yup = 2
ytics = 1      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Edit the following plot commands to correspond to data
# read in from data file

plot file1 using 1:2 with lines ls 1 title "dt"
replot file1 using 1:5 with lines ls 4 title "1/Rmax"
replot file1 using 1:3 with lines ls 9 title "1/Rmin"

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_stiffstuff.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

#*******************************************************************************************************************************************************#
#*************************************************************** ASY CHECK *****************************************************************************#
# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log kdt' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": mass fraction"
set title ds textcolor rgb title_color

file8 = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -5
yup = 5
ytics = 1.0      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Edit the following plot commands to correspond to data
# read in from data file

plot file8 using 1:10 with lines ls 2 title "4He"
replot file8 using 1:12 with lines ls 3 title "12C"
replot file8 using 1:14 with lines ls 1 title "16O"
replot file8 using 1:16 with lines ls 4 title "20Ne"
replot file8 using 1:18 with lines ls 5 title "24Mg"
replot file8 using 1:20 with lines ls 6 title "28Si"
replot file8 using 1:22 with lines ls 7 title "32S"
replot file8 using 1:24 with lines ls 8 title "36Ar"
replot file8 using 1:26 with lines ls 9 title "40Ca"
replot file8 using 1:28 with lines ls 10 title "44Ti"
replot file8 using 1:30 with lines ls 11 title "48Cr"
replot file8 using 1:32 with lines ls 12 title "52Fe"
replot file8 using 1:34 with lines ls 13 title "56Ni"
replot file8 using 1:36 with lines ls 14 title "60Zn"
replot file8 using 1:38 with lines ls 1 title "64Ge"
replot file8 using 1:40 with lines ls 2 title "68Se"


# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log kdt' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_ASY.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

#*****************************************************************************************************************************************************************#

#*************************************************************** PLOT DIFFC ******************************************************************************#

set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log diffC' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="ASY+PE T9=5 rho=1e8"
ds = ds.": diffC"
set title ds textcolor rgb title_color

fileC = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -16
yup = 0
ytics = 1.0      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Edit the following plot commands to correspond to data
# read in from data file

plot fileC using 1:8 with lines ls 1 title "|sumX - 1.0|"

#-----------------------------------------------------
# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log diffC' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_dC.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

quit

