# Sample 2D Vector Field 

# Plot (Irrotational) Vector Field
# See http://gnuplot.10905.n7.nabble.com/Vector-Fields-td3627.html
# See http://lavica.fesb.hr/cgi-bin/info2html?(gnuplot)arrow for arrowheads

set key inside left top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 600, 600   # x and y resolution of plot; larger is better
#unset xtics
#unset ytics

set style line 1 lt rgb "#5ea2c6" lw 1  # Define linestyle 1; rgb is hex #RRGGBB
set style line 2 lt rgb "#867961" lw 1  # Define linestyle 2; rgb is hex #RRGGBB

# Plot control parameters
hr = 5    # ranges
ss = 31   # points to evaluate vectors
vl = 0.25   # normalization factor for length

# X-Y range for plot
set xrange [-hr:hr]
set yrange [-hr:hr]

# Integer x-cordinates
set samples ss

# Integer y-cordinates
set isosamples ss

# The vector plot for gnuplot was originally designed for plotting discrete data from a file.
# The special filename "++" can be used to make vector plots from analytical expressions
# without actually having to produce the files.

# If at coordinates (x,y) the vector has x-projection dx and y-projection dy
# the appropriate plot command is of the form 
#     plot using "++" using 1:2:(expression for dx):(expression for dy) with vectors
# where in the expressions for dx and dy one references the variable x by $1 and the variable
# y by $2. This acts as if it were plotting a four-column data file, with the values in the
# four columns given by the four quantities separated by colons in the above plot command

# Simplest plot (no normalizations)
#plot "++" using 1:2:(-$2/($1**2+$2**2)):($1/($1**2+$2**2)) with vectors

# Plot normalized to unit length and the normalized length scaled by a factor of vl (the backslash
# is a line continuation, since a gnuplot command must be on a single line). The size and nature of
# the arrowhead on the vectors is controlled by "head filled size 0.10,20", which specifies a filled
# arrowhead with a length 0.10 that is proportional to the size of the x-axis (by default; it can
# be overridden) and 20 specifies that the arrow edges make a 20-degree angle with the line
# defining the vector.  The length of the vector line is set by the scale factor vl defined above.

# Counterclockwise flow
plot "++" using 1:2:(-$2*sqrt($1**2+$2**2)/($1**2+$2**2)*vl):($1*sqrt($1**2+$2**2)/($1**2+$2**2)*vl)\
with vectors head filled size 0.10,20 ls 1

# Clockwise flow
#plot "++" using 1:2:($2*sqrt($1**2+$2**2)/($1**2+$2**2)*vl):(-$1*sqrt($1**2+$2**2)/($1**2+$2**2)*vl)\
#with vectors head filled size 0.10,20 ls 1

# Create a postscript file of the plot

# Width and height of postscript figure in inches
width = 8
height = 5

set out "vectorField2D.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 32
replot               # Plot to postscript file

# Plot to PNG file

set out "vectorField2D.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2 enhanced font 'Arial,28'
replot
