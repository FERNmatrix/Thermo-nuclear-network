set xlabel 'Re x' 
set ylabel 'Im x'
set zlabel 'f(x)'
set xrange [-1.5:1.5]    
set yrange [-1.5:1.5]    
set zrange [-6.0:-0.0]

set isosamples 300       # Number surface mesh lines; larger -> slower and larger .ps
#set isosamples 500      # Production quality mesh
set hidden3d             # Remove lines that would be hidden if surface opaque
set key outside          # Move the legend outside the plot
#set palette model RGB   # or CMY (RGB default); 'help palette' other options
# Following define Palette interpolating between user-defined colors
#set palette defined (0 '#ff0000', 1 '#dddddd')  
set palette defined (0 'black', 1 'gold')  
set pm3d                  # Shaded color surface
set view 20, 33, 1, 2.0   # Viewing angles and scale factors
#set contour base         # Contours on base, surface, or both
#set cntrparam levels 10  # Target number of contour levels

a = -10  # Negative for mexican hat
b = 9

# Make 3D plot

splot a*(x**2+y**2) + 0.5*b*(x**2+y**2)**2 notitle

# Create a postscript file of the plot

# Width and height of postscript figure in inches
width = 8
height = 5

set out "example3D.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 32
replot               # Plot to postscript file

# Plot to PNG file

set out "example3D.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2 enhanced font 'Arial,28'
replot

