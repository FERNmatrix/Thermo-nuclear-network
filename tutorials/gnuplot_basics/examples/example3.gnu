# Multiplot example

# NOTE: the usual conversion to postscript seems only to
# pick up the final figure of the three plotted here, so
# maybe not so useful.

set multiplot layout 1,3 title "Multiplot layout for 3 plots"
#set xtics rotate        # Rotates x tic numbers by 90 degrees
set bmargin 2
#
set title "Plot 1"
unset key                # Suppresses inset key box for curves
plot sin(x)
#
set title "Plot 2"
unset key
plot cos(x)
#
set title "Plot 3"
unset key
plot x**2
#
unset multiplot

