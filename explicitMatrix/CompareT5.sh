#!/bin/bash

#Go to gnu out folder
cd gnu_out/

#copy data to Comparisons folder and enter folder
cp *.data Comparisons/T9\=5/
cd Comparisons/T9\=5/

#enter in gnuplot
gnuplot
