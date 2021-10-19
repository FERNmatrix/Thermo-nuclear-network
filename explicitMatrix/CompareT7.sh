#!/bin/bash

#Go to gnu out folder
cd gnu_out/

#copy data to Comparisons folder and enter folder
cp gnufile.data Comparisons/T9\=7/
cd Comparisons/T9\=7/

#enter in gnuplot
gnuplot
