 # Note:  next 2 lines of set key are part of single line

set key inside left top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000

set samples 600, 600   # x and y resolution of plot; larger is better

set style line 1 lt rgb "#5ea2c6" lw 1  # Define linestyle 1; rgb is hex #RRGGBB
set style line 2 lt rgb "#bb6255" lw 3  # Define linestyle 2; rgb is hex #RRGGBB

# Define a function to plot
gl0 = 0.05
a = -1.0
b = 1.0
gl(x) = gl0 + a*x**2 + b*x**4
gl2(x) = -a*x**2

set xlabel 'Order Parameter' 
set ylabel 'Value'
set xrange [-1.2:1.2]
set yrange [-0.25:0.3]

plot gl(x) ls 1, gl2(x) ls 2   # Plot gl w/linestyle 1 and gl2 w/linestyle 2
