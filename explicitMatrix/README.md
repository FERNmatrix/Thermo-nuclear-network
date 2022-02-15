1. Note that in Linux may need to add -lm compile flag to get the math.h header to work
for defining powf, expf, logf, ...  For example,
    
    gcc explicitMatrix.cpp -o explicitMatrix -lgsl -lgslcblas -lm -lstdc++

See https://www.includehelp.com/c-programming-questions/error-undefined-reference-to-pow-in-linux.aspx

2.  Code compiled above can then be executed with Linux with

./explicitMatrix  | tee temp.txt

where | tee temp.txt allows output to screen and also piped to a file temp.txt.

Execution for Mac or PC will depend on the C compiler installed on your machine.


                                 ******************************CLEAN SLATE BRANCH*************************************************
This branch will contain all impovements/corrections to original explicitMatrix code in MASTER
There will be a clean version of the working code titled explcitMatrix.cpp and then a version or 2 of the code with diagnostics and Autoscript capabailities

Currently (2/15/22):
	EMATS has up to date version of Timestepper with all corrections to PE, ASY and renormalization of X's.
		Currently ASY works fine, PE has some errors and leads to discrepency btwn the X's.
		The diagnostic plots for ASY kdt and diffX are commented out, but can easily be used if necessary
			The GNUplots gnuscript automates the plotting process for 5 plots
			The gnuDiagnostics.gnu file is the gnuscript to be used if ASYkdt and DiffX are wanting to be plotted
			
	NEM is old and the timestep is out of date, but a good thing to do with that is copy the current EMATS file and include the
	ASYkdt and diffC diagnostics as well as any ohter diagnosic tools we had before that could be useful.
	
	Makefile contains auto-compilation of 2 codes: EMATS and NEM
	
	CompareT7.sh has the automated commands to copy data to the compare/T9=7 folder and open gnuplot, using command: 
	<load 'GNUplots.gnu'> will load 5 plots
		1. dt vs t
		2. sumX vs t
		3. frac equil vs t
		4. dt & limiting rates vs t
		5. Mass fractions vs t
	CompareT5.sh does the same thing but for T9=5

