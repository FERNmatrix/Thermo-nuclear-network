Attached is a new version explicitMatrix.cpp, to which I have begun to add GSL, GSL_BLAS, and  C++ features, so it needs to be compiled as C++ with appropriate flags to link to GSL and BLAS libraries.  On my Fedora linux system using the gcc compiler this requires the library linking flags 

gcc explicitMatrix.cpp -o explicitMatrix -lgsl -lgslcblas -lm -lstdc++
where

-lgsl links GSL libraries
-lgslcblas links GSL BLAS libraries
-lm links support for Math.h header (required on my system; maybe not yours?)
-lstdc++ invokes the C++ compiler

Since C++ is a superset of C, all of the original C code still works.  The program explicitMatrix.cpp reads the same data files as before, so you should put this version in a directory containing the original data files.

A rough overview of the new features;

1.  Sample demonstration of of GSL vectors, and a specific construction of reaction vectors for all reactions in the network (see Section 3.2.3 of the notes "Matrix Formulation of Algebraically-Stabilized Explicit Integration for Thermonuclear Networks" sent earlier, hereafter referred to as "the notes").  Mostly in lines 430-450, the function makeReactionVectors ( ) beginning in line 1418, and the function compareGSLvectors(gsl_vector*, gsl_vector*) beginning in line 1608.  The array of reaction vectors rv[ ] is accessed primarily by setting the pointer *rvPt equal to the base address for the array and incrementing the pointer to step through the elements of the array.

2.  Demonstration of organizing data corresponding to the information for each isotope (proton number, neutron number, abundance, ...) into a C struct called networkData, which is defined in lines 248-280 and 371-372, with a pointer *networkDataPtr to the struct defined in line 283.  Its usage is shown in lines 1110-1159 of the function readNetwork(*fileName), where I demonstrate how to use the pointer *networkDataPtr to read data directly from the input file into the struct networkData (and also how to copy it from the struct to arrays and to the C++ Species object defined below). Also lines 459-489 illustrate retrieving and printing info from the struct.

3.  Implementation of a C++ class Species in lines 286-356 and its instantiation as an array of Species objects isotope[], with each element of the array a Species object corresponding to a single species (isotope) of the network.  These objects contain the same data as the struct networkData, but since they are derived from a C++ class the class encapsulates both the data and functions that operate on the data.  Now I read the data into the struct as described above and then copy into the Species objects isotope[ ], but we can employ the same pointer techniques used to read into the struct networkData to read the data directly into the objects isotope[ ].  (Pointers for C++ classes behave essentially as for C structs, except that they can point to functions in addition to data fields for a class;  examples are shown in this code.)

I wanted to give a general demonstration of using a struct.  However, C structs are like C++ classes but with much less functionality, so I think that we should use systematically C++ classes instead of C structs to encapsulate data and functions.

4. A classification of all reactions of the network into reaction group classes (see Sections 3.2.3 - 3.2.8 in the notes), which is essential bookkeeping for implementing the partial equilibrium approximation.  This uses the GSL API to define reaction vectors in makeReactionVectors( ) and compares them using compareGSLvectors( ) in sortReactionGroups( ) to assign a reaction group class to each reaction in the network (two reactions are in the same reaction group if their reaction vectors are equal, or equal up to an overall sign).

5.  Lines 1634-1692 implement a few new utility functions that will be handy at various places in developing the code.  Comments indicate what they do.
