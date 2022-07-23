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
set style line 14 lc rgb 'dark-grey' pt 4   # open circle
set style line 15 lc rgb 'dark-cyan' pt 4   # open circle
set style line 16 lc rgb 'purple' pt 4   # open circle
set style line 17 lc rgb 'orchid' pt 4   # open circle
set style line 18 lc rgb 'brown' pt 4   # open circle
set style line 19 lc rgb 'dark-green' pt 4   # open circle
set style line 20 lc rgb 'magenta' pt 4   # open circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color #font "Arial,22"
set ylabel 'Log X' textcolor rgb tic_color #font "Arial,22"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key right top outside    # Place legend inside top  
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="C++ T9=7 Asy rho=1e8 (no pf)"
ds = ds.": mass fraction"
set title ds textcolor rgb title_color

file1 = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = -5
xtics = 1    # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -14
yup = 0
ytics = 1    # Space between major y ticmarks
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

# Main CNO cycle

#plot file1 using 1:8 with lines ls 2 title "1H"
#replot file1 using 1:9 with lines ls 3 title "4He"
#replot file1 using 1:10 with lines ls 1 title "12C"
#replot file1 using 1:11 with lines ls 10 title "13C"
#replot file1 using 1:12 with lines ls 5 title "13N"
#replot file1 using 1:13 with lines ls 6 title "14N"
#replot file1 using 1:14 with lines ls 15 title "15N"
#replot file1 using 1:15 with lines ls 4 title "15O"

# Full CNO cycle (bi-cycle)

#plot file1 using 1:8 with lines ls 2 title "n"
#replot file1 using 1:9 with lines ls 3 title "1H"
#replot file1 using 1:10 with lines ls 1 title "4He"
#replot file1 using 1:11 with lines ls 10 title "12C"
#replot file1 using 1:12 with lines ls 5 title "13C"
#replot file1 using 1:13 with lines ls 6 title "13N"
#replot file1 using 1:14 with lines ls 15 title "14N"
#replot file1 using 1:15 with lines ls 1 title "15N"
#replot file1 using 1:16 with lines ls 2 title "15O"
#replot file1 using 1:17 with lines ls 3 title "16O"
#replot file1 using 1:18 with lines ls 4 title "17O"
#replot file1 using 1:19 with lines ls 5 title "18O"
#replot file1 using 1:20 with lines ls 6 title "17F"
#replot file1 using 1:21 with lines ls 7 title "18F"
#replot file1 using 1:22 with lines ls 8 title "19F"
#replot file1 using 1:23 with lines ls 9 title "20Ne"

# pp chains

#plot file1 using 1:8 with lines ls 2 lw 1.5 title "1H"
#replot file1 using 1:9 with lines ls 3 lw 1.5 title "2H"
#replot file1 using 1:10 with lines ls 1 lw 1.5 title "3He"
#replot file1 using 1:11 with lines ls 10 lw 1.5 title "4He"
#replot file1 using 1:12 with lines ls 5 lw 1.5 title "7Li"
#replot file1 using 1:13 with lines ls 6 lw 1.5 title "7Be"
#replot file1 using 1:14 with lines ls 15 lw 1.5 title "8B"


# Alpha network (X)

plot file1 using 1:8 with lines ls 2 lw 1.0 dashtype 1 title "4He"
replot file1 using 1:9 with lines ls 3 lw 1.0 dashtype 1 title "12C"
replot file1 using 1:10 with lines ls 1 lw 1.0 dashtype 1 title "16O"
replot file1 using 1:11 with lines ls 4 lw 1.0 dashtype 1 title "20Ne"
replot file1 using 1:12 with lines ls 5 lw 1.0 dashtype 1 title "24Mg"
replot file1 using 1:13 with lines ls 6 lw 1.0 dashtype 1 title "28Si"
replot file1 using 1:14 with lines ls 7 lw 1.0 dashtype 1 title "32S"
replot file1 using 1:15 with lines ls 8 lw 1.0 dashtype 1 title "36Ar"
replot file1 using 1:16 with lines ls 9 lw 1.0 dashtype 1 title "40Ca"
replot file1 using 1:17 with lines ls 10 lw 1.0 dashtype 1 title "44Ti"
replot file1 using 1:18 with lines ls 11 lw 1.0 dashtype 1 title "48Cr"
replot file1 using 1:19 with lines ls 12 lw 1.0 dashtype 1 title "52Fe"
replot file1 using 1:20 with lines ls 13 lw 1.0 dashtype 1 title "56Ni"
replot file1 using 1:21 with lines ls 14 lw 1.0 dashtype 1 title "60Zn"
replot file1 using 1:22 with lines ls 1 lw 1.0 dashtype 1 title "64Ge"
replot file1 using 1:23 with lines ls 2 lw 1.0 dashtype 1 title "68Se"


# 48-isotope network

#plot file1 using 1:8 with lines ls 2 lw 1 dashtype 2 title "(0)"
#replot file1 using 1:9 with lines ls 3 lw 1 dashtype 2 title "(1)"
#replot file1 using 1:10 with lines ls 1 lw 1 dashtype 2 title "(2)"
#replot file1 using 1:11 with lines ls 4 lw 1 dashtype 2 title "(3)"
#replot file1 using 1:12 with lines ls 5 lw 1 dashtype 2 title "(4)"
#replot file1 using 1:13 with lines ls 6 lw 1 dashtype 2 title "(5)"
#replot file1 using 1:14 with lines ls 7 lw 1 dashtype 2 title "(6)"
#replot file1 using 1:15 with lines ls 8 lw 1 dashtype 2 title "(7)"
#replot file1 using 1:16 with lines ls 9 lw 1 dashtype 2 title "(8)"
#replot file1 using 1:17 with lines ls 10 lw 1 dashtype 2 title "(9)"

#replot file1 using 1:18 with lines ls 11 lw 1 dashtype 2 title "(10)"
#replot file1 using 1:19 with lines ls 12 lw 1 dashtype 2 title "(11)"
#replot file1 using 1:20 with lines ls 13 lw 1 dashtype 2 title "(12)"
#replot file1 using 1:21 with lines ls 14 lw 1 dashtype 2 title "(13)"
#replot file1 using 1:22 with lines ls 3 lw 1 dashtype 2 title "(14)"
#replot file1 using 1:23 with lines ls 4 lw 1 dashtype 2 title "(15)"
#replot file1 using 1:24 with lines ls 15 lw 1 dashtype 2 title "(16)"
#replot file1 using 1:25 with lines ls 16 lw 1 dashtype 2 title "(17)"
#replot file1 using 1:26 with lines ls 17 lw 1 dashtype 2 title "(18)"
#replot file1 using 1:27 with lines ls 18 lw 1 dashtype 2 title "(19)"

#replot file1 using 1:28 with lines ls 19 lw 1 dashtype 2 title "(20)"
#replot file1 using 1:29 with lines ls 20 lw 1 dashtype 2 title "(21)"
#replot file1 using 1:30 with lines ls 5 lw 1 dashtype 2 title "(22)"
#replot file1 using 1:31 with lines ls 6 lw 1 dashtype 2 title "(23)"
#replot file1 using 1:32 with lines ls 7 lw 1 dashtype 2 title "(24)"
#replot file1 using 1:33 with lines ls 8 lw 1 dashtype 2 title "(25)"
#replot file1 using 1:34 with lines ls 9 lw 1 dashtype 2 title "(26)"
#replot file1 using 1:35 with lines ls 10 lw 1 dashtype 2 title "(27)"
#replot file1 using 1:36 with lines ls 11 lw 1 dashtype 2 title "(28)"
#replot file1 using 1:37 with lines ls 12 lw 1 dashtype 2 title "(29)"

#replot file1 using 1:38 with lines ls 13 lw 1 dashtype 2 title "(30)"
#replot file1 using 1:39 with lines ls 14 lw 1 dashtype 2 title "(31)"
#replot file1 using 1:39 with lines ls 13 lw 1 dashtype 2 title "(32)"
#replot file1 using 1:40 with lines ls 14 lw 1 dashtype 2 title "(33)"
#replot file1 using 1:41 with lines ls 15 lw 1 dashtype 2 title "(34)"
#replot file1 using 1:42 with lines ls 16 lw 1 dashtype 2 title "(35)"
#replot file1 using 1:43 with lines ls 17 lw 1 dashtype 2 title "(36)"
#replot file1 using 1:44 with lines ls 18 lw 1 dashtype 2 title "(37)"
#replot file1 using 1:45 with lines ls 19 lw 1 dashtype 2 title "(38)"
#replot file1 using 1:46 with lines ls 20 lw 1 dashtype 2 title "(39)"

#replot file1 using 1:47 with lines ls 21 lw 1 dashtype 2 title "(40)"
#replot file1 using 1:48 with lines ls 22 lw 1 dashtype 2 title "(41)"
#replot file1 using 1:49 with lines ls 1 lw 1 dashtype 2 title "(42)"
#replot file1 using 1:50 with lines ls 2 lw 1 dashtype 2 title "(43)"
#replot file1 using 1:51 with lines ls 3 lw 1 dashtype 2 title "(44)"
#replot file1 using 1:52 with lines ls 4 lw 1 dashtype 2 title "(45)"
#replot file1 using 1:53 with lines ls 5 lw 1 dashtype 2 title "(46)"
#replot file1 using 1:54 with lines ls 6 lw 1 dashtype 2 title "(47)"


# 70-isotope network

#plot file1 using 1:8 with lines ls 2 lw 1 dashtype 2 title "(0)"
#replot file1 using 1:9 with lines ls 3 lw 1 dashtype 2 title "(1)"
#replot file1 using 1:10 with lines ls 1 lw 1 dashtype 2 title "(2)"
#replot file1 using 1:11 with lines ls 4 lw 1 dashtype 2 title "(3)"
#replot file1 using 1:12 with lines ls 5 lw 1 dashtype 2 title "(4)"
#replot file1 using 1:13 with lines ls 6 lw 1 dashtype 2 title "(5)"
#replot file1 using 1:14 with lines ls 7 lw 1 dashtype 2 title "(6)"
#replot file1 using 1:15 with lines ls 8 lw 1 dashtype 2 title "(7)"
#replot file1 using 1:16 with lines ls 9 lw 1 dashtype 2 title "(8)"
#replot file1 using 1:17 with lines ls 10 lw 1 dashtype 2 title "(9)"

#replot file1 using 1:18 with lines ls 11 lw 1 dashtype 2 title "(10)"
#replot file1 using 1:19 with lines ls 12 lw 1 dashtype 2 title "(11)"
#replot file1 using 1:20 with lines ls 13 lw 1 dashtype 2 title "(12)"
#replot file1 using 1:21 with lines ls 14 lw 1 dashtype 2 title "(13)"
#replot file1 using 1:22 with lines ls 3 lw 1 dashtype 2 title "(14)"
#replot file1 using 1:23 with lines ls 4 lw 1 dashtype 2 title "(15)"
#replot file1 using 1:24 with lines ls 15 lw 1 dashtype 2 title "(16)"
#replot file1 using 1:25 with lines ls 16 lw 1 dashtype 2 title "(17)"
#replot file1 using 1:26 with lines ls 17 lw 1 dashtype 2 title "(18)"
#replot file1 using 1:27 with lines ls 18 lw 1 dashtype 2 title "(19)"

#replot file1 using 1:28 with lines ls 19 lw 1 dashtype 2 title "(20)"
#replot file1 using 1:29 with lines ls 20 lw 1 dashtype 2 title "(21)"
#replot file1 using 1:30 with lines ls 5 lw 1 dashtype 2 title "(22)"
#replot file1 using 1:31 with lines ls 6 lw 1 dashtype 2 title "(23)"
#replot file1 using 1:32 with lines ls 7 lw 1 dashtype 2 title "(24)"
#replot file1 using 1:33 with lines ls 8 lw 1 dashtype 2 title "(25)"
#replot file1 using 1:34 with lines ls 9 lw 1 dashtype 2 title "(26)"
#replot file1 using 1:35 with lines ls 10 lw 1 dashtype 2 title "(27)"
#replot file1 using 1:36 with lines ls 11 lw 1 dashtype 2 title "(28)"
#replot file1 using 1:37 with lines ls 12 lw 1 dashtype 2 title "(29)"

#replot file1 using 1:38 with lines ls 13 lw 1 dashtype 2 title "(30)"
#replot file1 using 1:39 with lines ls 14 lw 1 dashtype 2 title "(31)"
#replot file1 using 1:39 with lines ls 13 lw 1 dashtype 2 title "(32)"
#replot file1 using 1:40 with lines ls 14 lw 1 dashtype 2 title "(33)"
#replot file1 using 1:41 with lines ls 15 lw 1 dashtype 2 title "(34)"
#replot file1 using 1:42 with lines ls 16 lw 1 dashtype 2 title "(35)"
#replot file1 using 1:43 with lines ls 17 lw 1 dashtype 2 title "(36)"
#replot file1 using 1:44 with lines ls 18 lw 1 dashtype 2 title "(37)"
#replot file1 using 1:45 with lines ls 19 lw 1 dashtype 2 title "(38)"
#replot file1 using 1:46 with lines ls 20 lw 1 dashtype 2 title "(39)"

#replot file1 using 1:47 with lines ls 21 lw 1 dashtype 2 title "(40)"
#replot file1 using 1:48 with lines ls 22 lw 1 dashtype 2 title "(41)"
#replot file1 using 1:49 with lines ls 1 lw 1 dashtype 2 title "(42)"
#replot file1 using 1:50 with lines ls 2 lw 1 dashtype 2 title "(43)"
#replot file1 using 1:51 with lines ls 3 lw 1 dashtype 2 title "(44)"
#replot file1 using 1:52 with lines ls 4 lw 1 dashtype 2 title "(45)"
#replot file1 using 1:53 with lines ls 5 lw 1 dashtype 2 title "(46)"
#replot file1 using 1:54 with lines ls 6 lw 1 dashtype 2 title "(47)"
#replot file1 using 1:55 with lines ls 6 lw 1 dashtype 2 title "(48)"
#replot file1 using 1:56 with lines ls 6 lw 1 dashtype 2 title "(49)"

#replot file1 using 1:57 with lines ls 13 lw 1 dashtype 2 title "(50)"
#replot file1 using 1:58 with lines ls 14 lw 1 dashtype 2 title "(51)"
#replot file1 using 1:59 with lines ls 13 lw 1 dashtype 2 title "(52)"
#replot file1 using 1:60 with lines ls 14 lw 1 dashtype 2 title "(53)"
#replot file1 using 1:61 with lines ls 15 lw 1 dashtype 2 title "(54)"
#replot file1 using 1:62 with lines ls 16 lw 1 dashtype 2 title "(55)"
#replot file1 using 1:63 with lines ls 17 lw 1 dashtype 2 title "(56)"
#replot file1 using 1:64 with lines ls 18 lw 1 dashtype 2 title "(57)"
#replot file1 using 1:65 with lines ls 19 lw 1 dashtype 2 title "(58)"
#replot file1 using 1:66 with lines ls 20 lw 1 dashtype 2 title "(59)"

#replot file1 using 1:67 with lines ls 21 lw 1 dashtype 2 title "(60)"
#replot file1 using 1:68 with lines ls 22 lw 1 dashtype 2 title "(61)"
#replot file1 using 1:69 with lines ls 1 lw 1 dashtype 2 title "(62)"
#replot file1 using 1:70 with lines ls 2 lw 1 dashtype 2 title "(63)"
#replot file1 using 1:71 with lines ls 3 lw 1 dashtype 2 title "(64)"
#replot file1 using 1:72 with lines ls 4 lw 1 dashtype 2 title "(65)"
#replot file1 using 1:73 with lines ls 5 lw 1 dashtype 2 title "(66)"
#replot file1 using 1:74 with lines ls 6 lw 1 dashtype 2 title "(67)"
#replot file1 using 1:75 with lines ls 6 lw 1 dashtype 2 title "(68)"
#replot file1 using 1:76 with lines ls 6 lw 1 dashtype 2 title "(69)"



# Select isotopes from 150 network

#plot file1 using 1:8 with lines ls 2 lw 3 title "4He"
#replot file1 using 1:9 with lines ls 3 title "12C"
#replot file1 using 1:10 with lines ls 1 title "16O"
#replot file1 using 1:11 with lines ls 4 title "20Ne"
#replot file1 using 1:12 with lines ls 5 title "24Mg"
#replot file1 using 1:13 with lines ls 6 lw 3 title "28Si"
#replot file1 using 1:14 with lines ls 7 title "32S"
#replot file1 using 1:15 with lines ls 8 title "36Ar"
#replot file1 using 1:16 with lines ls 9 title "40Ca"
#replot file1 using 1:17 with lines ls 10 title "44Ti"
#replot file1 using 1:18 with lines ls 11 title "48Cr"
#replot file1 using 1:19 with lines ls 12 lw 3 title "52Fe"
#replot file1 using 1:20 with lines ls 13 lw 3 title "56Ni"
#replot file1 using 1:21 with lines ls 14 title "60Zn"

#replot file1 using 1:22 with lines ls 3 lw 3 title "1H"
#replot file1 using 1:23 with lines ls 4 lw 3 title "n"
#replot file1 using 1:24 with lines ls 15 title "13C"
#replot file1 using 1:25 with lines ls 16 title "14N"
#replot file1 using 1:26 with lines ls 17 title "29Si"
#replot file1 using 1:27 with lines ls 18 title "31P"
#replot file1 using 1:28 with lines ls 19 title "64Zn"
#replot file1 using 1:29 with lines ls 20 title "60Ni"
#replot file1 using 1:30 with lines ls 5 title "55Co"
#replot file1 using 1:31 with lines ls 6 title "24Al"

#replot file1 using 1:32 with lines ls 7 title "19Fl"
#replot file1 using 1:33 with lines ls 8 title "22Na"
#replot file1 using 1:34 with lines ls 9 title "22Ne"
#replot file1 using 1:35 with lines ls 10 title "23Mg"
#replot file1 using 1:36 with lines ls 11 title "7Be"
#replot file1 using 1:37 with lines ls 12 lw 3 title "14O"
#replot file1 using 1:38 with lines ls 13 lw 3 title "17O"
#replot file1 using 1:39 with lines ls 14 title "24Al"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,18"
set key top right font "Arial,14"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,21"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,21"

# Plot to postscript file

set out "gnuplot_X.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 18
replot               # Plot to postscript file

# Plot to PNG file

set out "gnuplot_X.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2
replot

quit
