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
width = 12
height = 12

# x-axis resolution
set samples 1000

# Line styles.  
# For lines: plot x with lines
# For points: plot x with points

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
set style line 11 lc rgb 'red' pt 6   # open circle
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

set pointsize 1.0    # Size of the plotted points

set key right top outside font "Arial,8"    # Place legend inside top 
unset key            # Don't show legend

set timestamp       # Date/time

ds="C++ Asy+PE 70 iso, 200 plot steps"
ds = ds.": T9=6 rho=1e8"
set title noenhanced   # Symbols like underscore not interpreted as markup
set title ds textcolor rgb title_color


# -------- Axis ranges and ticmarks -----------

xlow = -12
xup = -5
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -14
yup = 0
ytics = 2    # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------------------

file1 = "plot1.data"

# Reference data

#refFile =  "dataRef/gnufile_150_viktorProfile_400_asyRef_c++.data"
#refFile = "dataRef/gnufile_alpha_T9_5_1e7_asy.data"
#refFile = "dataRef/gnufile_alpha_T9_7_1e8_asy_C++_PF.data"
#refFile = "dataRef/gnufile_alpha_victorProfile_400_asyRef_c++.data"  
#refFile = "dataRef/nova125D_sumX_1.000.data"  
#refFile = "dataRef/gnuplot_alpha_viktorProfileSmooth_asyRef_c++.data" 
#refFile = "dataRef/gnuplot_alpha_viktorProfileSmooth_asyRef_java.data"
refFile = "dataRef/gnufile_70_T9=6_rho=1e8_asyRef_c++.data"


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

#plot file1 using 1:8 with lines ls 11 lw 1.0 dashtype 1 title "4He"
#replot file1 using 1:9 with lines ls 3 lw 1.0 dashtype 1 title "12C"
#replot file1 using 1:10 with lines ls 1 lw 1.0 dashtype 1 title "16O"
#replot file1 using 1:11 with lines ls 4 lw 1.0 dashtype 1 title "20Ne"
#replot file1 using 1:12 with lines ls 5 lw 1.0 dashtype 1 title "24Mg"
#replot file1 using 1:13 with lines ls 6 lw 1.0 dashtype 1 title "28Si"
#replot file1 using 1:14 with lines ls 7 lw 1.0 dashtype 1 title "32S"
#replot file1 using 1:15 with lines ls 8 lw 1.0 dashtype 1 title "36Ar"
#replot file1 using 1:16 with lines ls 9 lw 1.0 dashtype 1 title "40Ca"
#replot file1 using 1:17 with lines ls 10 lw 1.0 dashtype 1 title "44Ti"
#replot file1 using 1:18 with lines ls 11 lw 1.0 dashtype 1 title "48Cr"
#replot file1 using 1:19 with lines ls 12 lw 1.0 dashtype 1 title "52Fe"
#replot file1 using 1:20 with lines ls 13 lw 1.0 dashtype 1 title "56Ni"
#replot file1 using 1:21 with lines ls 14 lw 1.0 dashtype 1 title "60Zn"
#replot file1 using 1:22 with lines ls 1 lw 1.0 dashtype 1 title "64Ge"
#replot file1 using 1:23 with lines ls 2 lw 1.0 dashtype 1 title "68Se"


# Select isotopes from the 70-isotope network

plot file1 using 1:10 with lines ls 0 lw 2 dashtype 1 title "(2 4He)"
replot file1 using 1:12 with lines ls 1 lw 2 dashtype 1 title "(4 12C)"
replot file1 using 1:16 with lines ls 2 lw 2 dashtype 1 title "(8 160)"
replot file1 using 1:20 with lines ls 4 lw 2 dashtype 1 title "(12 20Ne)"
replot file1 using 1:25 with lines ls 5 lw 2 dashtype 1 title "(17 Mg24)"
replot file1 using 1:32 with lines ls 6 lw 2 dashtype 1 title "(24 28Si)"
replot file1 using 1:41 with lines ls 7 lw 2 dashtype 1 title "(33 32S)"
replot file1 using 1:49 with lines ls 8 lw 2 dashtype 1 title "(41 36Ar)"
replot file1 using 1:54 with lines ls 9 lw 2 dashtype 1 title "(46 40Ca)"
replot file1 using 1:57 with lines ls 10 lw 2 dashtype 1 title "(49 44Ti)"

replot file1 using 1:61 with lines ls 11 lw 2 dashtype 1 title "(53 48Cr)"
replot file1 using 1:70 with lines ls 12 lw 2 dashtype 1 title "(62 54Fe)"
replot file1 using 1:71 with lines ls 13 lw 2 dashtype 1 title "(63 55Fe )"
replot file1 using 1:74 with lines ls 14 lw 2 dashtype 1 title "(66 57Co)"
replot file1 using 1:77 with lines ls 3 lw 2 dashtype 1 title "(69 58Ni)"
#replot file1 using 1:23 with lines ls 4 lw 2 dashtype 1 title "(15)"
#replot file1 using 1:24 with lines ls 15 lw 2 dashtype 1 title "(16)"
#replot file1 using 1:25 with lines ls 16 lw 2 dashtype 1 title "(17)"
#replot file1 using 1:26 with lines ls 17 lw 2 dashtype 1 title "(18)"
#replot file1 using 1:27 with lines ls 18 lw 2 dashtype 1 title "(19)"

# Reference data, select isotopes from the 70-isotope network

replot refFile using 1:10 with lines ls 0 lw 2 dashtype 2 title "(2 4He ref)"
replot refFile using 1:12 with lines ls 1 lw 2 dashtype 2 title "(4 12C ref)"
replot refFile using 1:16 with lines ls 2 lw 2 dashtype 2 title "(8 160 ref)"
replot refFile using 1:20 with lines ls 4 lw 2 dashtype 2 title "(12 20Ne ref)"
replot refFile using 1:25 with lines ls 5 lw 2 dashtype 2 title "(17 Mg24 ref)"
replot refFile using 1:32 with lines ls 6 lw 2 dashtype 2 title "(24 28Si ref)"
replot refFile using 1:41 with lines ls 7 lw 2 dashtype 2 title "(33 32S ref)"
replot refFile using 1:49 with lines ls 8 lw 2 dashtype 2 title "(41 36Ar ref)"
replot refFile using 1:54 with lines ls 9 lw 2 dashtype 2 title "(46 40Ca ref)"
replot refFile using 1:57 with lines ls 10 lw 2 dashtype 2 title "(49 44Ti ref)"

replot refFile using 1:61 with lines ls 11 lw 2 dashtype 2 title "(53 48Cr ref)"
replot refFile using 1:70 with lines ls 12 lw 2 dashtype 2 title "(62 54Fe ref)"
replot refFile using 1:71 with lines ls 13 lw 2 dashtype 2 title "(63 55Fe  ref)"
replot refFile using 1:74 with lines ls 14 lw 2 dashtype 2 title "(66 57Co ref)"
replot refFile using 1:77 with lines ls 3 lw 2 dashtype 2 title "(69 58Ni ref)"
#replot file1 using 1:23 with lines ls 4 lw 2 dashtype 2 title "(15 ref)"
#replot file1 using 1:24 with lines ls 15 lw 2 dashtype 2 title "(16 ref)"
#replot file1 using 1:25 with lines ls 16 lw 2 dashtype 2 title "(17 ref)"
#replot file1 using 1:26 with lines ls 17 lw 2 dashtype 2 title "(18 ref)"
#replot file1 using 1:27 with lines ls 18 lw 2 dashtype 2 title "(19 ref)"



# Following lines output generated in plotfile1 --> plot1.data by
#
#
#    for(int i=0; i<maxPlotIsotopes; i++){
#        plotXlist[i] = i;
#    }
#
# in plotFileSetup() of explicitMatrix.cpp.  Following will plot
# all isotopes for maxPlotIsotopes <= 150 set in explicitMatrix.cpp.
# Numbers in parentheses are the isotope index of the species in the 
# calculation, so this can plot generically for any network, but the 
# meaning of the numbers varies from network to network. The species
# (isotope) number for each isotope is output to the screen at the 
# end of a explicitMatrix.cpp calculation.

#plot file1 using 1:8 with lines ls 0 lw 1 dashtype 1 title "(0)"
#replot file1 using 1:9 with lines ls 1 lw 1 dashtype 1 title "(1)"
#replot file1 using 1:10 with lines ls 2 lw 1 dashtype 1 title "(2)"
#replot file1 using 1:11 with lines ls 4 lw 1 dashtype 1 title "(3)"
#replot file1 using 1:12 with lines ls 5 lw 1 dashtype 1 title "(4)"
#replot file1 using 1:13 with lines ls 6 lw 1 dashtype 1 title "(5)"
#replot file1 using 1:14 with lines ls 7 lw 1 dashtype 1 title "(6)"
#replot file1 using 1:15 with lines ls 8 lw 1 dashtype 1 title "(7)"
#replot file1 using 1:16 with lines ls 9 lw 1 dashtype 1 title "(8)"
#replot file1 using 1:17 with lines ls 10 lw 1 dashtype 1 title "(9)"

#replot file1 using 1:18 with lines ls 11 lw 1 dashtype 1 title "(10)"
#replot file1 using 1:19 with lines ls 12 lw 1 dashtype 1 title "(11)"
#replot file1 using 1:20 with lines ls 13 lw 1 dashtype 1 title "(12)"
#replot file1 using 1:21 with lines ls 14 lw 1 dashtype 1 title "(13)"
#replot file1 using 1:22 with lines ls 3 lw 1 dashtype 1 title "(14)"
#replot file1 using 1:23 with lines ls 4 lw 1 dashtype 1 title "(15)"
#replot file1 using 1:24 with lines ls 15 lw 1 dashtype 1 title "(16)"
#replot file1 using 1:25 with lines ls 16 lw 1 dashtype 1 title "(17)"
#replot file1 using 1:26 with lines ls 17 lw 1 dashtype 1 title "(18)"
#replot file1 using 1:27 with lines ls 18 lw 1 dashtype 1 title "(19)"

#replot file1 using 1:28 with lines ls 19 lw 1 dashtype 1 title "(20)"
#replot file1 using 1:29 with lines ls 20 lw 1 dashtype 1 title "(21)"
#replot file1 using 1:30 with lines ls 5 lw 1 dashtype 1 title "(22)"
#replot file1 using 1:31 with lines ls 6 lw 1 dashtype 1 title "(23)"
#replot file1 using 1:32 with lines ls 7 lw 1 dashtype 1 title "(24)"
#replot file1 using 1:33 with lines ls 8 lw 1 dashtype 1 title "(25)"
#replot file1 using 1:34 with lines ls 9 lw 1 dashtype 1 title "(26)"
#replot file1 using 1:35 with lines ls 10 lw 1 dashtype 1 title "(27)"
#replot file1 using 1:36 with lines ls 11 lw 1 dashtype 1 title "(28)"
#replot file1 using 1:37 with lines ls 12 lw 1 dashtype 1 title "(29)"

#replot file1 using 1:38 with lines ls 13 lw 1 dashtype 1 title "(30)"
#replot file1 using 1:39 with lines ls 14 lw 1 dashtype 1 title "(31)"
#replot file1 using 1:39 with lines ls 13 lw 1 dashtype 1 title "(32)"
#replot file1 using 1:40 with lines ls 14 lw 1 dashtype 1 title "(33)"
#replot file1 using 1:41 with lines ls 15 lw 1 dashtype 1 title "(34)"
#replot file1 using 1:42 with lines ls 16 lw 1 dashtype 1 title "(35)"
#replot file1 using 1:43 with lines ls 17 lw 1 dashtype 1 title "(36)"
#replot file1 using 1:44 with lines ls 18 lw 1 dashtype 1 title "(37)"
#replot file1 using 1:45 with lines ls 19 lw 1 dashtype 1 title "(38)"
#replot file1 using 1:46 with lines ls 20 lw 1 dashtype 1 title "(39)"

#replot file1 using 1:47 with lines ls 21 lw 1 dashtype 1 title "(40)"
#replot file1 using 1:48 with lines ls 22 lw 1 dashtype 1 title "(41)"
#replot file1 using 1:49 with lines ls 1 lw 1 dashtype 1 title "(42)"
#replot file1 using 1:50 with lines ls 2 lw 1 dashtype 1 title "(43)"
#replot file1 using 1:51 with lines ls 3 lw 1 dashtype 1 title "(44)"
#replot file1 using 1:52 with lines ls 4 lw 1 dashtype 1 title "(45)"
#replot file1 using 1:53 with lines ls 5 lw 1 dashtype 1 title "(46)"
#replot file1 using 1:54 with lines ls 6 lw 1 dashtype 1 title "(47)"
#replot file1 using 1:55 with lines ls 6 lw 1 dashtype 1 title "(48)"
#replot file1 using 1:56 with lines ls 6 lw 1 dashtype 1 title "(49)"

#replot file1 using 1:57 with lines ls 13 lw 1 dashtype 1 title "(50)"
#replot file1 using 1:58 with lines ls 14 lw 1 dashtype 1 title "(51)"
#replot file1 using 1:59 with lines ls 13 lw 1 dashtype 1 title "(52)"
#replot file1 using 1:60 with lines ls 14 lw 1 dashtype 1 title "(53)"
#replot file1 using 1:61 with lines ls 15 lw 1 dashtype 1 title "(54)"
#replot file1 using 1:62 with lines ls 16 lw 1 dashtype 1 title "(55)"
#replot file1 using 1:63 with lines ls 17 lw 1 dashtype 1 title "(56)"
#replot file1 using 1:64 with lines ls 18 lw 1 dashtype 1 title "(57)"
#replot file1 using 1:65 with lines ls 19 lw 1 dashtype 1 title "(58)"
#replot file1 using 1:66 with lines ls 20 lw 1 dashtype 1 title "(59)"

#replot file1 using 1:67 with lines ls 21 lw 1 dashtype 1 title "(60)"
#replot file1 using 1:68 with lines ls 22 lw 1 dashtype 1 title "(61)"
#replot file1 using 1:69 with lines ls 1 lw 1 dashtype 1 title "(62)"
#replot file1 using 1:70 with lines ls 2 lw 1 dashtype 1 title "(63)"
#replot file1 using 1:71 with lines ls 3 lw 1 dashtype 1 title "(64)"
#replot file1 using 1:72 with lines ls 4 lw 1 dashtype 1 title "(65)"
#replot file1 using 1:73 with lines ls 5 lw 1 dashtype 1 title "(66)"
#replot file1 using 1:74 with lines ls 6 lw 1 dashtype 1 title "(67)"
#replot file1 using 1:75 with lines ls 6 lw 1 dashtype 1 title "(68)"
#replot file1 using 1:76 with lines ls 6 lw 1 dashtype 1 title "(69)"

#replot file1 using 1:77 with lines ls 2 lw 1 dashtype 1 title "(70)"
#replot file1 using 1:78 with lines ls 3 lw 1 dashtype 1 title "(71)"
#replot file1 using 1:79 with lines ls 1 lw 1 dashtype 1 title "(72)"
#replot file1 using 1:80 with lines ls 4 lw 1 dashtype 1 title "(73)"
#replot file1 using 1:81 with lines ls 5 lw 1 dashtype 1 title "(74)"
#replot file1 using 1:82 with lines ls 6 lw 1 dashtype 1 title "(75)"
#replot file1 using 1:83 with lines ls 7 lw 1 dashtype 1 title "(76)"
#replot file1 using 1:84 with lines ls 8 lw 1 dashtype 1 title "(77)"
#replot file1 using 1:85 with lines ls 9 lw 1 dashtype 1 title "(78)"
#replot file1 using 1:86 with lines ls 10 lw 1 dashtype 1 title "(79)"

#replot file1 using 1:87 with lines ls 11 lw 1 dashtype 1 title "(80)"
#replot file1 using 1:88 with lines ls 12 lw 1 dashtype 1 title "(81)"
#replot file1 using 1:89 with lines ls 13 lw 1 dashtype 1 title "(82)"
#replot file1 using 1:90 with lines ls 14 lw 1 dashtype 1 title "(83)"
#replot file1 using 1:91 with lines ls 3 lw 1 dashtype 1 title "(84)"
#replot file1 using 1:92 with lines ls 4 lw 1 dashtype 1 title "(85)"
#replot file1 using 1:93 with lines ls 15 lw 1 dashtype 1 title "(86)"
#replot file1 using 1:94 with lines ls 16 lw 1 dashtype 1 title "(87)"
#replot file1 using 1:95 with lines ls 17 lw 1 dashtype 1 title "(88)"
#replot file1 using 1:96 with lines ls 18 lw 1 dashtype 1 title "(89)"

#replot file1 using 1:97 with lines ls 19 lw 1 dashtype 1 title "(90)"
#replot file1 using 1:98 with lines ls 20 lw 1 dashtype 1 title "(91)"
#replot file1 using 1:99 with lines ls 5 lw 1 dashtype 1 title "(92)"
#replot file1 using 1:100 with lines ls 6 lw 1 dashtype 1 title "(93)"
#replot file1 using 1:101 with lines ls 7 lw 1 dashtype 1 title "(94)"
#replot file1 using 1:102 with lines ls 8 lw 1 dashtype 1 title "(95)"
#replot file1 using 1:103 with lines ls 9 lw 1 dashtype 1 title "(96)"
#replot file1 using 1:104 with lines ls 10 lw 1 dashtype 1 title "(97)"
#replot file1 using 1:105 with lines ls 11 lw 1 dashtype 1 title "(98)"
#replot file1 using 1:106 with lines ls 12 lw 1 dashtype 1 title "(99)"

#replot file1 using 1:107 with lines ls 13 lw 1 dashtype 1 title "(100)"
#replot file1 using 1:108 with lines ls 14 lw 1 dashtype 1 title "(101)"
#replot file1 using 1:109 with lines ls 13 lw 1 dashtype 1 title "(102)"
#replot file1 using 1:110 with lines ls 14 lw 1 dashtype 1 title "(103)"
#replot file1 using 1:111 with lines ls 15 lw 1 dashtype 1 title "(104)"
#replot file1 using 1:112 with lines ls 16 lw 1 dashtype 1 title "(105)"
#replot file1 using 1:113 with lines ls 17 lw 1 dashtype 1 title "(106)"
#replot file1 using 1:114 with lines ls 18 lw 1 dashtype 1 title "(107)"
#replot file1 using 1:115 with lines ls 19 lw 1 dashtype 1 title "(108)"
#replot file1 using 1:116 with lines ls 20 lw 1 dashtype 1 title "(109)"

#replot file1 using 1:117 with lines ls 21 lw 1 dashtype 1 title "(110)"
#replot file1 using 1:118 with lines ls 22 lw 1 dashtype 1 title "(111)"
#replot file1 using 1:119 with lines ls 1 lw 1 dashtype 1 title "(112)"
#replot file1 using 1:120 with lines ls 2 lw 1 dashtype 1 title "(113)"
#replot file1 using 1:121 with lines ls 3 lw 1 dashtype 1 title "(114)"
#replot file1 using 1:122 with lines ls 4 lw 1 dashtype 1 title "(115)"
#replot file1 using 1:123 with lines ls 5 lw 1 dashtype 1 title "(116)"
#replot file1 using 1:124 with lines ls 6 lw 1 dashtype 1 title "(117)"
#replot file1 using 1:125 with lines ls 6 lw 1 dashtype 1 title "(118)"
#replot file1 using 1:126 with lines ls 6 lw 1 dashtype 1 title "(119)"

#replot file1 using 1:127 with lines ls 11 lw 1 dashtype 1 title "(120)"
#replot file1 using 1:128 with lines ls 12 lw 1 dashtype 1 title "(121)"
#replot file1 using 1:129 with lines ls 13 lw 1 dashtype 1 title "(122)"
#replot file1 using 1:130 with lines ls 14 lw 1 dashtype 1 title "(123)"
#replot file1 using 1:131 with lines ls 3 lw 1 dashtype 1 title "(124)"
#replot file1 using 1:132 with lines ls 4 lw 1 dashtype 1 title "(125)"
#replot file1 using 1:133 with lines ls 15 lw 1 dashtype 1 title "(126)"
#replot file1 using 1:134 with lines ls 16 lw 1 dashtype 1 title "(127)"
#replot file1 using 1:135 with lines ls 17 lw 1 dashtype 1 title "(128)"
#replot file1 using 1:136 with lines ls 18 lw 1 dashtype 1 title "(129)"

#replot file1 using 1:137 with lines ls 13 lw 1 dashtype 1 title "(130)"
#replot file1 using 1:138 with lines ls 14 lw 1 dashtype 1 title "(131)"
#replot file1 using 1:139 with lines ls 13 lw 1 dashtype 1 title "(132)"
#replot file1 using 1:140 with lines ls 14 lw 1 dashtype 1 title "(133)"
#replot file1 using 1:141 with lines ls 15 lw 1 dashtype 1 title "(134)"
#replot file1 using 1:142 with lines ls 16 lw 1 dashtype 1 title "(135)"
#replot file1 using 1:143 with lines ls 17 lw 1 dashtype 1 title "(136)"
#replot file1 using 1:144 with lines ls 18 lw 1 dashtype 1 title "(137)"
#replot file1 using 1:145 with lines ls 19 lw 1 dashtype 1 title "(138)"
#replot file1 using 1:146 with lines ls 20 lw 1 dashtype 1 title "(139)"

#replot file1 using 1:147 with lines ls 21 lw 1 dashtype 1 title "(140)"
#replot file1 using 1:148 with lines ls 22 lw 1 dashtype 1 title "(141)"
#replot file1 using 1:149 with lines ls 1 lw 1 dashtype 1 title "(142)"
#replot file1 using 1:150 with lines ls 2 lw 1 dashtype 1 title "(143)"
#replot file1 using 1:151 with lines ls 3 lw 1 dashtype 1 title "(144)"
#replot file1 using 1:152 with lines ls 4 lw 1 dashtype 1 title "(145)"
#replot file1 using 1:153 with lines ls 5 lw 1 dashtype 1 title "(146)"
#replot file1 using 1:154 with lines ls 6 lw 1 dashtype 1 title "(147)"
#replot file1 using 1:155 with lines ls 6 lw 1 dashtype 1 title "(148)"
#replot file1 using 1:156 with lines ls 6 lw 1 dashtype 1 title "(149)"


# Reset font sizes for .eps and .png output2

set timestamp font "Arial,16"

set title ds textcolor rgb title_color font "Arial,18"
set key top right font "Arial,14"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,22"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,22"

# Plot to postscript file

set out "gnuplot_X.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Symbol,22"
replot               # Plot to postscript file

# Plot to PNG file

#set out "gnuplot_X.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit
