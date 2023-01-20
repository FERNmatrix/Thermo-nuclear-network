/*
 * Code explicitMatrix.cpp to implement explicit algebraic integration of astrophysical 
 * thermonuclear networks. Execution assuming use of Fedora Linux and GCC compiler: Compile with
 * 
 *     gcc explicitMatrix.cpp -o explicitMatrix -O2 -lgsl -lgslcblas -lm -lstdc++
 * 
 * [may need to install gsl-devel development package (on Fedora this required
 * dnf install gsl-devel-2.4-8.fc30.x86_64) if it can't find gsl headers in the compile.]
 * In this compile command the -o specifies the name of the executable created in the
 * compile, the -lgsl flag links to GSL libraries, the -O2 flag sets the level of optimiztion
 * for the compiler, the -lgslcblas flag
 * links to GSL BLAS libraries, the -lm flag may be required in GCC Linux to get the 
 * math.h header to work for defining pow, exp, log,... (see https://
 * www.includehelp.com/c-programming-questions/error-undefined-reference-to-pow-in-linux.aspx),
 * and lstdc++ is the link flag specifying the C++ compiler to use. 
 * Alternatively you can use the makefile defined the directory to compile with
 * 
 *      make
 * 
 * at the command line. Note that the present make file uses the g++ rather than gcc compiler.
 * 
 * If you plan to use the GDB debugger, add a -g flag:
 * 
 *     gcc explicitMatrix.cpp -o explicitMatrix -O2 -lgsl -lgslcblas -lm -lstdc++ -g
 * 
 * The compiled code created above can be executed with
 * 
 *     ./explicitMatrix | tee temp.txt
 * 
 * where | tee temp.txt is unix shell script outputting to screen and also piping to a file temp.txt. 
 * If using GDB debug mode, set the -g flag as above for the compiler and execute in debug mode with
 * 
 *      gdb explicitMatrix | tee debug.out
 * 
 * Further information about using GDB may be found in the directory DEBUGGER.
 * 
 * Memory checks with valgrind:
 * 
 * valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./explicitMatrix 2>&1 | tee valgrind_out.txt
 * 
 * Compile explicitMatrix.cpp with the -g flag to get valgrind output with human readable line numbers. Further information about valgrind may be found in the directory VALGRIND.
 * 
 * Execution for other Linux systems, or Mac or PC, will depend on the C/C++ compiler installed on 
 * your machine but should be similar to above.
 *  
 * ----------------------------------
 * To set up a specific calculation:
 * ----------------------------------
 * 
 * 1. Change values of ISOTOPES and SIZE.
 * 2. Change input files for networkFile and rateLibraryFile.
 * 3. Change doASY, doQSS, and doPE to choose Asy, Asy+PE, QSS, QSS+PE options.
 * 4. Change control parameters like stop_time, massTol,...
 * 5. Change values of T9_start and rho_start if constant T and rho (hydroProfile=false).
 * 6. If using hydro profile: hydroProfile=true, set HydroFile[], set maxHydroEntries,
 *    choose true or false for plotHydroProfile
 * ----------------------------------
 * 
 * 
 * AUTHORS:
 * ---------------
 * Nick Brey
 * Ashton DeRousse
 * Adam Cole
 * Raghav Chari
 * Jay Billings
 * Mike Guidry
 * ----------------
 */


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;
//using std::string;

// Define some CPU timing utilities. Usage:
//
//     START_CPU;
//     ... code to be timed ...
//     STOP_CPU;
//     PRINT_CPU
//
// These may be used to time CPU processes.

clock_t startCPU,stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("Timer: %7.4e sec used",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU (fprintf(pFile,"in %7.4e seconds\n",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU2 (fprintf(pFile2,"in %7.4e seconds\n",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU3 (fprintf(plotfile1,"in %7.4e seconds\n",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU4 (fprintf(plotfile2,"in %7.4e seconds\n",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPUD (fprintf(pFileD,"in %g seconds\n",(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define PRINT_CPU_TEST (printf("\nTimer Test: %g ms used by CPU\n",1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));  

// Function signatures:

void devcheck(int);
void readLibraryParams(char *);
void readNetwork(char *);
void readhydroProfile(char *);
void writeNetwork(void);
void writeRates(char *);
void writeAbundances(void);
void setRG(int,int,int);
void setSpeciesfplus(int,double);
void setSpeciesfminus(int,double);
void setSpeciesdYdt(int,double);
void assignRG(void);
void getmaxdYdt(void);
void restoreEquilibriumProg(void);
void evolveToEquilibrium(void);
bool isoIsInRG(int,int);
void computeReactionFluxes();
void updateTheFluxes();
void updateY0(void);
void showY(void);
void showParameters(void);
double dE_halfstep(void);
void sumFplusFminus(void);
void restoreBe8(void);
void plotFileSetup(void);
void toPlotNow(void);
void updatePF(void);
void displayRGstatus(void);
double returnNetworkMass(void);
void networkMassDifference(void);
bool checkForCNOCycle(void);
void correctCNOCycle(void);
void indexCNOCycle(void);

/*
 * ------------------------------------------------------------------------------------------
 * SOME SAMPLE NETWORKS:
 * 
 * Network   ISOTOPES     SIZE     networkFile[]               rateLibraryFile[]
 * 3=alpha          3        8     data/network_3alpha.inp     data/rateLibrary_3alpha.data
 * 4-alpha          4       14     data/network_4alpha.inp     data/rateLibrary_4alpha.data
 * alpha           16       48     data/network_alpha.inp      data/rateLibrary_alpha.data
 * pp               7       28     data/network_pp.inp         data/rateLibrary_pp.data
 * main cno         8       22     data/network_cno.inp        data/rateLibrary_cno.data
 * full cno        16      134     data/network_cnoAll.inp     data/rateLibrary_cnoAll.data
 * 48              48      299     data/network_48.inp         data/rateLibrary_48.data
 * 70(C-O)         70      598     data/network_70.inp         data/rateLibrary_70.data
 * 70(4He)         70      598     data/network_70_alpha.inp   data/rateLibrary_70.data
 * 116            116     1135     data/network_116.inp        data/rateLibrary_116.data
 * nova134        134     1566     data/network_nova134.inp    data/rateLibrary_nova134.data
 * 150 (12C-16O)  150     1604     data/network_150.inp        data/rateLibrary_150.data
 * 150 (solar)    150     1604     data/network_150_solar.inp  data/rateLibrary_150.data
 * 150 (he4)      150     1604     data/network_150_he4.inp    data/rateLibrary_150.data
 * 194            194     2232     network_194.inp             data/rateLibrary_194.data
 * 268            268     3175     network_268.inp             data/rateLibrary_268.data
 * 365 (12C-16O)  365     4395     data/network_365.inp        data/rateLibrary_365.data
 * 365 (solar)    365     4395     data/network_365_solar.inp  data/rateLibrary_365.data
 * tidalSN_alpha   16       48     data/network_alpha_he4.inp  data/rateLibrary_alpha.data
 * big bang         8       64     data/network_bigbang.inp    data/rateLibrary_bigbang.data
 * 28              28      104     data/network_28.inp         data/rateLibrary_28.data
 * 30P             47      283     data/network_test30P.inp    data/rateLibrary_test30P.data
 * ------------------------------------------------------------------------------------------
 */


//  SIZE defines the number of reactions to be calculated. ISOTOPES defines the number 
//  of isotopes in each network.  These sizes are hardwired for now but eventually we may want 
//  to read them in and assign them dynamically.

#define ISOTOPES 134                   // Max isotopes in network (e.g. 16 for alpha network)
#define SIZE 1566                      // Max number of reactions (e.g. 48 for alpha network)

#define plotSteps 200                // Number of plot output steps
#define LABELSIZE 35                  // Max size of reaction string a+b>c in characters
#define PF 24                         // Number entries partition function table for isotopes
#define THIRD 0.333333333333333
#define TWOTHIRD 0.66666666666667
#define ECON 9.5768e17                // Convert MeV/nucleon/s to erg/g/s
#define LOG10 0.434294481903251       // Conversion natural log to log10
#define MEV 931.494                   // Conversion of amu to MeV
#define GZ 1.0e-24                    // Constant to ensure 1/max(num, GZ) never divides by 0

#define unitd static_cast<double>(1.0)  // Constant double equal to 1
#define zerod static_cast<double>(0.0)  // Constant double equal to 0

// File pointers for diagnostics output. Corresponding filenames declared 
// at top of main.

FILE* pFileD;
FILE* pfnet;

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format 
// for this file. Standard filenames for test cases are listed in table above.

char networkFile[] = "data/network_nova134.inp";

// Filename for input rates library data. The file rateLibrary.data output by 
// the Java code through the stream toRateData has the expected format for this 
// file.  Standard filenames for test cases are listed in table above.

char rateLibraryFile[] = "data/rateLibrary_nova134.data";

// Whether to use constant T and rho (hydroProfile false), in which case a
// constant T9 = T9_start and rho = rho_start are used, or to read
// in a hydrodynamical profile of T and rho versus time (hydroProfile true),
// in which case the file to be read in is specified by the character variable 
// hydroFile[].

bool hydroProfile = true; 

// Filename for input file containing a hydro profile in temperature
// and density that is used if hydroProfile = true. Sample hydro profile 
// files included in the data subdirectory are
//
//    data/torch47Profile.data         // Very hot Type Ia supernova zone
//    data/nova125DProfile.inp         // Representative zone in nova explosion
//    data/tidalSNProfile_100.inp      // Zone in tidal supernova explosion
//
// Use SplineInterpolator to interpolate in table read in. If hydroProfile and 
// plotHydroProfile are true, the hydro profile used for the temperature and 
// density in the calculation is also output to the file gnu_out/hydroProfileInput.data
// in format suitable for gnuplot.

char hydroFile[] = "data/nova125DProfile_400.inp"; //"data/rosswog.profile";

// Control output of hydro profile (if one is used) to plot file.

static const bool plotHydroProfile = true;

const static int maxHydroEntries = 403; //2622; // Max entries hydro profile

// Control printout of flux data (true to print, false to suppress).
// Lots of data, so most useful for small networks.
 
static const bool plotFluxes = false;

// Plot output controls and file pointers

static const int maxPlotIsotopes = min(ISOTOPES, 365);   // # species to plot
int plotXlist[maxPlotIsotopes];           // Array of species plot indices

// Pointers to data output files

FILE* plotfile1;
FILE* plotfile2;
FILE* plotfile3;
FILE* plotfile4;
FILE* plotfile5;

// Control flags for diagnostic output to files. Note that setting showDetails
// or showDetails2 true may generate large output files (MB to GB for large networks).

bool showAddRemove = true;   // Show addition/removal of RG from equilibrium
bool showDetails = false;    // Controls diagnostics to pFileD -> gnu_out/diagnostics.data
bool showDetails2 = false;   // Controls diagnostics to pfnet -> gnu_out/network.data

// Control which explicit algebraic approximations are used. Eventually
// this should be set from a data file. To use asymptotic set doASY true
// (which toggles doQSS to false). To use quasi-steady-state (QSS),set 
// doASY false (which toggles doQSS to true). doPE can be true or false 
// with either Asymptotic or QSS. The boolean showPE allows display of the number
// of reaction groups (RG) that would be in equilibrium if PE approximation were
// being implemented. It is true only if Asy or QSS, but PE not being
// implemented.

bool doASY = true;            // Whether to use asymptotic approximation
bool doQSS = !doASY;          // Whether to use QSS approximation 
bool doPE = false;             // Implement partial equilibrium also
bool showPE = !doPE;          // Show RG that would be in equil if doPE=false

string intMethod = "";        // String holding integration method
string ts;                    // Utility string

// Temperature and density variables. Temperature and density can be
// either constant, or read from a hydro profile as a function of time.

double T9;                    // Current temperature in units of 10^9 K
double rho;                   // Current density in units of g/cm^3

// EfromMasses = false if energy plotted from Q values; EfromMasses = true if energy 
// is plotted from masses by weighing the network before and after timestep.
// Generally should be set to EfromMasses = true, since energy from Q values and
// mass differences are essentially the same for no PE approximation, but
// for PE approximation the fluxes set to zero are not quite zero in reality
// so fluxes computed from Q values are increasingly in error as more and
// more reaction groups come into equilibrium.  Conversely, energies 
// computed by weighing the network (mass differences weighted by abundances)
// are approximately in agreement with exact calculations at the same level
// as the agreement of the abundances. For now we will compute the energies 
// both ways but by setting EfromMasses = true we will ensure that the more
// correct energies (E and dE/dt) computed from mass differences are output to 
// the plotting streams.

bool EfromMasses = true;    // Set to true except for testing purposes       

// Energy variables (from Q values)

double ERelease;              // Total energy released
double dERelease;             // Energy released per unit time
double netdERelease;          // Energy released in timestep

// Energy variables (from mass differences)

double networkMass;           // Total mass in network at end of timestep
double lastNetworkMass;       // Total network mass at beginning of timestep
double netdEReleaseA;         // Net energy release in timestep dE from mass diffs
double dEReleaseA;            // Differential E release dE/dt over timestep from masses
double EReleaseA;             // Total E released up to this time from masses

// Partition function controls. If dopf = true, reaction rates are
// corrected by temperature-dependent partition functions.  However
// partition function factors differ from 1 only at high temperature
// so we only implement partition function corrections if T9 > pfCut9,
// where pfCut9 is a cutoff temperature in units of T9. Typically in
// realistic calculation we would choose dopf = true and pfCut9 = 1.0.

bool dopf = true;
double pfCut9 = 1.0;

// Temperatures in units of 10^9 K for partition function table (see pf[]
// in the class Species). 

double Tpf[PF] = {0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    
// Array holding the value of the partition function for each isotope at
// the current temperature if dopf = true and T9 > pfCut9.
    
double currentPF[ISOTOPES];

// Array to hold whether given species satisfies asymptotic condition
// True (1) if asyptotic; else false (0).

bool isAsy[ISOTOPES];         // True if isotope is asymptotic
double asycheck;              // Species asymptotic if asycheck > 1.0
double asyFrac = 0.0;         // Fraction isotopes that are asymptotic

// Whether isotope part of any RG in partial equilibrium this timestep

bool isotopeInEquil[ISOTOPES]; 

// isotopeInEquil[] from last timestep

bool isotopeInEquilLast[ISOTOPES]; 

// Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
// realistic calculation the temperature and density will be passed from the hydro 
// code in an operator-split coupling of this network to hydro. Here we hardwire
// constant values for testing purposes, or read in a temperature and density
// hydro profile if hydroProfile is true.

double T9_start = 0.025;      // Initial temperature in units of 10^9 K
double rho_start = 100;        // Initial density in g/cm^3

// Integration time data. The variables start_time and stop_time 
// define the range of integration (all time units in seconds),
// and dt_start sets the initial integration timestep. In an operator-split 
// coupling  start_time will be 0, stop_time will correspond to the length
// of the hydro timestep and dt_init will likely be something like the 
// last timestep of the previous network integration (for the preceding 
// hydro timestep). Here we hardwire them for testing purposes.
// The variable startplot_time allows the plotting interval output
// in gnu_out/gnufile.data to be a subset of the full15 integration interval. 
// Generally, startplot_time > start_time.  By default the stop time for
// plotting is the same as the stop time for integration, stop_time.

double start_time = 1e-20;             // Start time for integration
double logStart = log10(start_time);   // Base 10 log start time
double startplot_time =  5e-6;         // Start time for plot output
double stop_time = 1e7;                // Stop time for integration
double logStop = log10(stop_time);     // Base-10 log stop time5
double dt_start = 0.01*start_time;     // Initial value of integration dt
double dt_saved;                       // Full timestep used for this int step
double t_saved;                        // Start time this timestep (end t for last step)
double dt_half;                        // Half of full timestep
double dt_change;                      // Change in proposed dt from last timestep
double t_end;                          // End time for this timestep
double dt_new;                         // Variable used in computeNextTimeStep()
double dtmin;                          // Variable used in computeNextTimeStep()
double dt_desired;                     // dt desired if not prevented by plot timestep
double dt_ceiling = 0.1;               // Max timestep is dt_ceiling*t, for accuracy

double dt_FE = dt_start;               // Max stable forward Euler timestep
double dt_EA = dt_start;               // Max asymptotic timestep

int dtMode;                            // Dual dt stage (0=full, 1=1st half, 2=2nd half)

double massTol_asy = 1e-7;             // Tolerance param if no reactions equilibrated
double massTol_asyPE = 6e-6;           // Tolerance param if some reactions equilibrated
double massTol = massTol_asy;          // Timestep tolerance parameter for integration
double downbumper = 0.7;               // Asy dt decrease factor
double sf = 1e25;                      // dt_FE = sf/fastest rate
int maxit = 100;                       // Max asy dt iterations allowed for a step
int iterations;                        // # iterations in step to conserve particles 
int totalIterations;                   // Total number of iterations, all steps til now
int mostIterationsPerStep = 0;         // Most iterations in a timestep
int maxIterationStep;                  // Step where mostIterationsPerStep occurred
double maxIterationTime;               // Time where mostIterationsPerStep occurred
double Error_Observed;                 // Observed integration error
double Error_Desired;                  // Desired max local integration error
double E_R;                            // Ratio actual to desired error
double EpsA = massTol_asy;             // Absolute error tolerance
double EpsR = 2.0e-4;                  // Relative error tolerance (not presently used)

// Apply cycle stabilization (CS) to CNO if X[H]<startX_fixCNO and fixCNO=true.

bool fixCNO = false;                    // Whether to apply cycle stabilization (CS) to CNO 
double startX_fixCNO = 2e-4;           // Fraction hydrogen mass fraction to start CS

bool CNOinNetwork = false;             // Whether currently applying CS correction
bool fixingCNO_now = false;            // Whether CS being applied at this timestep
double startX_fixCNO_time;             // Time when startX_fixCNO reached for H mass fraction

// Isotopic index for CNO isotopes. Set in function indexCNOcycle()
// if fixCNO is true.

int index12C;
int index13N;
int index13C;
int index14N;
int index15O;
int index15N;
int index1H;
int index4He;

// Reaction index for relevant CNO main cycle reactions. Each reaction has two
// components, so make 2D arrays initialized to -1. Set in function indexCNOcycle()
// if fixCNO is true.

int index14N_pgamma[] = {-1, -1};
int index13C_pgamma[] = {-1, -1};
int index15N_pgamma[] = {-1, -1};

// equilTime is time to begin imposing partial equilibrium if doPE=true. Hardwired but 
// eventually should be determined by the program.  In the Java version this was sometimes
// needed because starting PE test too early could lead to bad results.  This is 
// probably an error in the Java version, since if operating properly nothing should
// be changed at a timestep if nothing satisfies PE condition.  Thus, we should not need
// this in a final version for stability, but it might still be useful since early in
// a calculation typically nothing satisfies PE, so checking for it is a waste of time.
// On the other hand, the check costs little computing time so to make the code more
// universal it may be best to check for equilibration from the beginning of the 
// calculation. 

double equilTime = start_time;    // Time to begin checking for PE
double equiTol = 0.015;           // Tolerance for checking whether Ys in RG in equil
double deviousMax = 0.2;         // Max allowed deviation from equil k ratio in timestep
bool useDevious = false;          // Use thisDevious (true) of equil pops (false) to set equil
bool useEquilY = true;            // Use equilibrium values of Y to impose PE

double thisDevious;               // Deviation of kratio from equil
double mostDevious = 0.0;         // Largest current deviation of kratio from equil
int mostDeviousIndex;             // Index of RG with mostDevious
int choice1;                      // Diagnostic variable for new timestepper
int choice2;                      // Diagnostic variable for new timestepper

double logTnow;                   // Log10 of current temp
double logRhoNow;                 // Log10 of current rho
double dt;                        // Current integration timestep
double dtLast;                    // Last timestep
double t;                         // Current time in integration
int totalTimeSteps;               // Number of integration timesteps taken
int totalTimeStepsZero;           // Timestep when plotting starts
int plotCounter;                  // Plot output counter
double logTimeSpacing;            // Constant spacing of plot points in log time
int totalAsy;                     // Total number of asymptotic isotopes

double sumX;                      // Sum of mass fractions X(i).  Should be 1.0.
double sumXeq;                    // Sum X species in at least 1 equilibrated RG
double sumXNeq;                   // Sum X species not in an equilibrated RG
double sumXtot;                   // Sum X all species
double XcorrFac;                  // Equil normalization factor for timestep
double sumXlast;                  // sumX from last timestep
double diffX;                     // sumX - sumXlast
double diffXzero;                 // diffX before computeTimeStep_EA (a,b) iteration 
double diffXfinal;                // diffX after computeTimeStep_EA (a,b) iteration
double sumXtrue;                  // True (un-renormalized) sumX

double Rate[SIZE];                // Reaction rate from Reactions::computeRate()
double Flux[SIZE];                // Flux from Reactions::computeFlux()
bool reacLibWarn;                 // Warning flag for T out of ReacLib bounds
double warnTime;                  // Time for 1st violation associated with reacLibWarn
double warnT;                     // Temperature at warnTime
double maxdYdt;                   // Maximum current dY/dt in network
int maxdYdtIndex;                 // Isotopic index of species with max dY/dt


// --- Species data in following arrays also contained in fields of class Species

int Z[ISOTOPES];                  // Array holding Z values for isotopes
int N[ISOTOPES];                  // Array holding N values for isotopes
int AA[ISOTOPES];                 // Array holding A values for isotopes
double Y[ISOTOPES];               // Array holding current abundances Y for isotopes
double Y0[ISOTOPES];              // Array holding abundances at beginning of timestep
double Ystore[ISOTOPES];          // Array storing current abundances for later restore
double YfullStep[ISOTOPES];       // Array holding full-step Y update
double X[ISOTOPES];               // Array holding mass fractions X for isotopes
double massExcess[ISOTOPES];      // Array holding mass excesses for isotopes
const static int isoLen = 5;      // Max character length for isoLabel[][]
char isoLabel[ISOTOPES][isoLen];  // Isotope labels (max 5 characters; e.g. 238pu)
double dYDt[ISOTOPES];            // Rate of change for Y


// --- Reaction data in following arrays also contained in fields of class Reaction

char reacLabel[SIZE][LABELSIZE]; // Reaction labels (e.g. he4+c12-->o16) 
int RGclass[SIZE];               // Reaction Group class (PE) for reaction (1-5)
int RGMemberIndex[SIZE];         // Member index within its reaction group
string RGstring[SIZE];           // Schematic RG; e.g. a <->b+c
int NumReactingSpecies[SIZE];    // Number of reactant isotopes for the reaction
int NumProducts[SIZE];           // Number of product isotopes for the reaction
int reacZ[SIZE][4];              // Holds Z for each reactant isotope
int reacN[SIZE][4];              // Holds N for each reactant isotope
int prodZ[SIZE][4];              // Holds Z for each product isotope
int prodN[SIZE][4];              // Holds N for each product isotope
int ReactantIndex[SIZE][4];      // Index of isotope vector for reactant isotope
int ProductIndex[SIZE][4];       // Index of isotope vector for product isotope
int isPEforward[SIZE];           // Whether labeled "forward" reaction in PE scheme


int numberSpecies;               // Actual # species in network (usually = ISOTOPES)
int numberReactions;             // Actual # reactions in network (usually = SIZE)

// Array with entries positive integer if a reaction increases the population of the isotope 
// (contributes to F+), negative integer if it decreases it (contributes to F-) and 0 if the 
// reaction does not change the population of the isotope. This array is populated 
// in the function ReactionVector::parseF(). It characterizes the structure 
// of the network and thus has to be calculated only once for a given network.

int reacMask[ISOTOPES][SIZE]; 

// Define an array RV of type std::vector that will hold reaction vectors for
// the system.

vector <int> RV[SIZE];

// Total number of F+ and F- terms in the network

int totalFplus = 0;
int totalFminus = 0;

// Arrays to hold time,temperature,and density in hydro profile

int hydroLines;                    // Number of hydro profile lines read in
double hydroTime[maxHydroEntries];
double hydroTemp[maxHydroEntries];
double hydroRho[maxHydroEntries];

// Arrays to hold non-zero fluxes in the network. Corresponding memory will be 
// allocated dynamically below with malloc.

double* Fplus;               // Dynamical 1D array for non-zero F+ (Dim totalFplus)
double* Fminus;              // Dynamical 1D array for non-zero F- (Dim totalFminus)

double FplusZero[ISOTOPES];  // Used in QSS approximation

double keff[ISOTOPES];       // Effective decay constant for given isotope
double keffZero[ISOTOPES];   // Used in QSS approximation

// Arrays to hold number species factors for F+ and F- arrays. For example,for 12C+12C ->
// 4He + 20Ne the species factor is 1 for 4He and 20Ne diff. equation terms but 2 for
// 12C diff equation term (since the reaction involves one 4He and one 20Ne,but two 12C.

double* FplusFac;     // Dynamically allocated 1D array for number species factor in F+ terms
double* FminusFac;    // Dynamically allocated 1D array for number species factor in F- terms

double* FplusSum;     // Sum of F+ for each isotope
double* FminusSum;    // Sum of F- for each isotope
double* dF;           // Net flux FplusSum-FminusSum for each isotope

int* FplusMax;        // Upper index for each isotope in the Fplus array
int* FplusMin;        // Lower index for each isotope in the Fplus array
int* FminusMax;       // Upper index for each isotope in the Fminus array
int* FminusMin;       // Lower index for each isotope in the Fminus array

// Arrays to hold index of upper boundary for entries of each isotope in the Fplus and
// Fminus arrays. Corresponding memory will be allocated dynamically below with malloc

int* FplusIsotopeCut;    // Upper index for each isotope in Fplus (Dim numberSpecies)
int* FminusIsotopeCut;   // Upper index for each isotope in Fminus (Dim numberSpecies)

int* numFluxPlus;        // Number of finite F+ components for each isotope (Dim numberSpecies)
int* numFluxMinus;       // Number of finite F- components for each isotope (Dim numberSpecies)

int* FplusIsotopeIndex;  // Array containing the isotope index for each F+ (Dim totalFplus)
int* FminusIsotopeIndex; // Array containing the isotope index for each F- (Dim totalFminus)

/*
 * Arrays to hold the mapping of the Fplus and Fminus arrays of fluxes to the
 * master flux array Flux[SIZE]. For example,MapFplus[0] will hold the value of
 * the index i in Flux[i] that corresponds to the reaction that generates Fplus[0].
 * These will be allocated dynamically below.  In general a given reaction flux appears
 * in more than one Fplus or Fminus array.  So we compute the flux for each reaction only
 * once for a given temperature and density,and then use this mapping to assign a given
 * flux to its occurrence in multiple Fplus and Fminus entries without having to
 * recompute fluxes.
 */

int* MapFplus;    // Index mapper for Fplus (Dim totalFplus)
int* MapFminus;   // Index mapper for Fminus (Dim totalFminus)

// Arrays for temporary storage. Will be allocated dynamically below

int* tempInt1;
int* tempInt2;

char dasher[] = "---------------------------------------------";


// ---------------------------------
// Partial equilibrium quantities
// ---------------------------------

int numberRG;                 // Number of partial equilibrium reaction groups
int RGnumberMembers[SIZE];    // # members each RG; determined in class ReactionVectors
int numberSingletRG;          // # RG with only a single member reactions
double eqFrac = 0.0;          // Fraction of RG equilibrated
bool showRGstatus = false;     // Output RG status to pfnet -> gnu_out/network.data" 

// Define array to hold the ReactionGroup object index for each reaction. There 
// are n reaction groups in a network and each reaction belongs to one and only
// one reaction group.  RGindex[m] is the index (0,1,... n) of the reaction 
// group to which reaction m belongs. This array is populated by the function
// ReactionVector::sortReactionGroups(). The ReactionGroup object index is also
// contained in the rgindex field of the Reaction objects reaction[i]. 

int RGindex[SIZE];

// Could save a little memory by assigning dynamically to size [numberRG][ISOTOPES] 
// once numberRG determined

bool RGisoMembers[SIZE][ISOTOPES];

bool reacIsActive[SIZE];    // False if reaction has been removed by PE

int totalEquilReactions;    // Total equilibrated reactions for isotope
int totalEquilRG;           // Total equilibrated reaction groups
bool normPEX = true;        // Normalize X after PE evol (normally true)

gsl_matrix* fluxes;
gsl_vector* abundances;

// Temporary utility quantities to hold fastest and slowest rate data at
// a given timestep. The index quantities are the index of the
// corresponding reaction in reactions[].

double fastestCurrentRate = 0.0;
int fastestCurrentRateIndex;
double slowestCurrentRate = 1e30;
int slowestCurrentRateIndex;
double fastestOverallRate = 0.0;
int fastestOverallRateIndex;
double slowestOverallRate = 1e30;
int slowestOverallRateIndex;
double timeMaxRate = 0.0;

// Logic that will be used to convert all 8-Be two alpha particles 
// since lifetime for 8-Be to decay to two alpha particles is short 
// compared with typical integration timesteps.

bool hasBe8 = false;     // Whether network contains BE8
int indexBe8= -1;        // Species index of Be8 if hasBe8 = true
bool hasAlpha = false;   // Whether network contains He4
int indexAlpha = -1;     // Species index of He4 if hasAlpha = true


// -------------------------------------------------------------------------
// Arrays to hold output quantities at plot timesteps. The variable
// plotSteps is hardwired above,but eventually should be input
// at time of execution from data files.
// -------------------------------------------------------------------------

// Target plot times for plot steps. Computed from value of plotSteps in
// Utilities::log10Spacing().

double plotTimeTargets[plotSteps];           // Target plot times
double nextPlotTime;                         // Next plot output time



//----------------CLASS DEFINITIONS ----------------

/*
 Class SplineInterpolator to implement 1D cubic spline interpolation. 
 Adapted from algorithms in NumericaaddToEquilibrium()l Recipes. For 1D interpolations,use 
 the method spline to set up an interpolation table and then use the 
 method splint to interpolate in the independent variable.
 The class also makes available the utility function bisection,which 
 finds the indices in a 1-D table that bracket a particular entry in 
 the table (assuming the entries to increase monotonically). 
 
 REFERENCE:  
 
 https://en.wikipedia.org/wiki/Spline_interpolation
 http://fourier.eng.hmc.edu/e176/lectures/ch7/node6.html
 http://www.foo.be/docs-free/Numerical_Recipe_In_C/c3-3.pdf
 
 */


class SplineInterpolator{
    
private:
    
    int numberPoints; 
    double x[maxHydroEntries];
    double y[maxHydroEntries];
    double y2[maxHydroEntries];
    double u[maxHydroEntries];
    
public:
    
    // Constructor creates a SplineInterpolator object for the arrays
    // xarray and yarray passed using the pointers *xarray and *yarray.
    
    SplineInterpolator(int points, double* xarray, double* yarray) { 

        numberPoints = points;
        
    };
    
    
    /*------------------------------------------------------------------------------
     *   Adaptation of 1D cubic spline algorithm from Numerical Recipes.
     *   This method processes the array to compute and store the second derivatives
     *   that will be used to do later spline interpolation.  It assumes
     *   the existence of an array xarray of independent variables and an array 
     *   yarray(x) of dependent variables,with the xarray array monotonically 
     *   increasing.  The second derivatives are stored in the array y2.  
     * ------------------------------------------------------------------------------*/
    
    
    void spline(double* xarray, double* yarray, int size1, int size2) {
        
        int n = size1;  
        int m = size2; 

        if(size1 != size2){
            printf("\nWarning: lengths of x and y(x) arrays not equal:");
            printf("\nxarray_length=%d yarray_length=%d",size1,size2);
            printf("\nEXIT\n");
            exit(-1);
        }
                    
        // Copy arrays passed in constructor to internal arrays for later use.
        
        for(int i=0; i<size1; i++){
            x[i] = xarray[i];
            y[i] = yarray[i];
        }
        
        // Natural spline boundary conditions
            
        y2[0] = 0.0;
        u[0] = 0.0;
        double qn = 0.0;
        double un = 0.0;
        double signum;
        double sigden;
        double sig;
        double p;
        
        // Decomposition loop of tridiagonal algorithm
        
        for (int i=1; i<= n-2; i++) {
            signum = x[i] - x[i-1];
            sigden = x[i+1] - x[i-1];
            sig = signum/sigden;
            p = sig*y2[i-1] +2.0;
            y2[i] = (sig-1.0)/p;
            u[i] = ( 6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) 
                /(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1] )/p;
        }
        
        y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
        
        // Backsubstitution loop of tridiagonal algorithm
        
        for (int i = n-2; i >= 0; i--){
            y2[i] = y2[i]*y2[i+1] + u[i];
        }
    }
    
    /*-------------------------------------------------------------------------
     * Method splint calculates the 1D cubic spline polynomial for an arbitrary
     * argument xvalue once the second derivative table has been constructed
     * using the method spline.  This method returns the interpolated value 
     * y = f(xvalue) and may be called any number of times once the second
     * derivative table has been created.
     * ---------------------------------------------------------------------------*/
    
    double splint(double xvalue) {
        
        int n = numberPoints;
        
        // Return -1 with error message if argument out of table bounds
        
        if (xvalue < x[0] || xvalue > x[n-1]) {
            printf("\nArgument ( %7.4e ) Out of Table Bounds %7.4e to %7.4e",
                xvalue,x[0],x[n-1]);
            return -1;
        }
        
        // Call bisection method to bracket entry xvalue with indices ilow and ihigh
        
        int ilow = bisection(x,xvalue);      
        int ihigh = ilow + 1;                  
        double h = x[ihigh]-x[ilow];
        
        // Evaluate cubic spline polynomial and return interpolated value
        
        double a = (x[ihigh]-xvalue)/h;
        double b = (xvalue-x[ilow])/h;
        return a*y[ilow] + b*y[ihigh] 
            + ((a*a*a-a)*y2[ilow]+(b*b*b-b)*y2[ihigh])*h*h/6.0;
        
    }
    
    
    /*-----------------------------------------------------------------------------
     * For an array xarray[] and argument xvalue,public method bisection finds 
     * the indices ilow and ihigh=ilow+1 that bracket the position of xvalue 
     * in the array.  (Method returns the value of ilow,from which 
     * ihigh = ilow+1).  The array is assumed to be monotonically increasing
     * in value. Sample method call:
     * 
     *	double [] x = {10,12,14,16,18};
     *	double xvalue = 12.2;
     *	int ilow = instance.bisection(x,xvalue);
     *	int ihigh= ilow + 1;
     * 
     * In this case,the method returns ilow = 1 and therefore ihigh = 2.  The method
     * checks first that xvalue is within the bounds of the table and returns -1 if
     * this condition is not satisfied.
     *------------------------------------------------------------------------------*/
    
    
   int bisection(double* xarray, double xvalue){
        
        int n = numberPoints;
        
        // Check that xvalue is within bounds of the table.  If not,quit
        // with error message and return -1
        
        double minx = xarray[0];
        double maxx = xarray[n-1];
        if(xvalue > maxx || xvalue < minx){
            printf("Abort bisection: argument (%7.4e) Out of Bounds",xvalue);
            return -1;
        }
        
        int ilow = 0;
        int ihigh = n-1;
        int i;
        
        while( (ihigh-ilow > 1) ){
            i = (ihigh+ilow)/2;
            if(xarray[i] > xvalue){
                ihigh=i;
            } else {
                ilow=i;
            }
        }
        
        // ilow and ilow+1 now bracket xvalue
        
        return ilow;
    }
    
};




/* Class Utilities to hold utility useful utility functions.  Functions are
 * declared static so that they can be invoked without having to instantiate
 * objects of type Utilities.  For example,Utilities::returnNetIndexZN (Z,N).
 * We will often have other classes inherit from Utilities so that they can
 * access its static functions.  This class definition must precede definitions of 
 * classes that inherit from it.
 */


class Utilities{
    
    private:
    
    public:
        
        static const int MAXCSIZE = 4500;   // Max string size in stringToChar(string)
        
        // Static function Utilities::showTime() to return date and local time as a 
        // character array
        
        static char* showTime(){
            time_t now = time(0);         // Current date/time
            char* tnow = ctime(&now);     // convert to string
            return tnow;
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::log10Spacing() to find num equal log10 
        // spacings between two numbers start and stop. The equally log-spaced
        // numbers are placed in the array v passed with the pointer * v.
        // -------------------------------------------------------------------------
        
        static void log10Spacing(double start, double stop, int num, double* v){
            
            double logtmin = log10(start);
            double logtmax = log10(stop);
            double tempsum = logtmin;
            logTimeSpacing = (logtmax - logtmin) / (double) num;
            v[0] = pow(10, logtmin);
            
            double tlow = 0;
            double tup = 0;
            double dtmax = 0;
            
            for(int i=0; i<num; i++){

                tempsum += logTimeSpacing;
                v[i] = pow(10, tempsum);       // v[i] mapped to plotTimeTarget[]

            }
        }
        
        
        // Static function Utilities::plotHydroprofile() to send hydro profile 
        // to plotting file. Only invoked if hydroProfile and plotHydroProfile 
        // are true. This is the hydro profile read in, which is output to
        // the file gnu_out/hydroProfileInput.data. The values of temperature and
        // density interpolated during the calculation are output to the file
        // plot4.data.
        
        static void outputHydroProfile(){
            
            FILE* pHydro;
            pHydro = fopen("gnu_out/hydroProfileInput.data","w");
            if( pHydro == NULL ) {
                fprintf(stderr,"Couldn't open file: %s\n",strerror(errno));
                exit(1);
            }
            
            fprintf(pHydro,"#  time          T          rho");
            
            for (int i=0; i<hydroLines; i++){
                fprintf(pHydro,"\n%6.4e  %6.4e  %6.4e",
                    hydroTime[i],hydroTemp[i],hydroRho[i]);
            }
            
            fclose (pHydro);
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::sumMassFractions() to return the current sum of the
        // mass fractions X(i) in the network. If the network conserves particle
        // number this sum should be equal to 1.0.
        // -------------------------------------------------------------------------
        
        static double sumMassFractions(void) {
            
            double sum = 0.0;
            for(int i=0; i<ISOTOPES; i++){
                sum += X[i];
            }
            return sum; 
        }
        
        
        // ------------------------------------------------------------------
        // Static function Utilities::sumXEquil() to return the sum of mass 
        // fractions for isotopes participating in at least one RG currently
        // in partial equilibrium
        // ------------------------------------------------------------------
        
        static double sumXEquil() {
            
            double sum = 0;
            for(int i=0; i<ISOTOPES; i++){
                if( isotopeInEquil[i] ){
                    sum += X[i];
                }
            }
            return sum;
        }
        
        
        // -----------------------------------------------------------------------
        // Static function Utilities::sumXNotEquil() to return the sum of mass 
        // fractions for isotopes not participating in any RG currently in
        // partial equilibrium
        // -----------------------------------------------------------------------
        
        static double sumXNotEquil() {
            
            double sum = 0;
            for(int i=0; i<ISOTOPES; i++){
                if( !isotopeInEquil[i] ){
                    sum += X[i];
                }
            }
            return sum;
        }
        
    
        // -------------------------------------------------------------------------
        // Static function Utilities::returnNetIndexZN(Z,N) to return the network 
        // vector index given Z and N for the isotope.  Returns -1 if no match.
        // -------------------------------------------------------------------------
        
        static int returnNetIndexZN (int z, int n) {
            
            for (int i = 0; i < numberSpecies; i++) {
                if (Z[i] == z && N[i] == n) return i;
            }
            return -1;
        }
        
        
        // -----------------------------------------------------------------------
        // Static function Utilities::returnNetIndexSymbol(char* symbol) to 
        // return the network vector index given the symbol for the isotope.
        // -----------------------------------------------------------------------
        
        static int returnNetIndexSymbol(char* symbol) {
            
            int result;
            for (int i = 0; i < numberSpecies; i++) {
                result = strcmp(isoLabel[i], symbol);
                if (result == 0){
                    return i;
                }
            }
            return -1;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::isInNet(Z, N) to return true if given (Z, N) 
        // is in the network, false otherwise.
        // ----------------------------------------------------------------------
        
        static bool isInNet(int Z, int N) {
            
            if (returnNetIndexZN(Z, N) < 0) {
                return false;
            } else {
                return true;
            }
        }
        
        // ----------------------------------------------------------------------
        // Static function Utilities::returnReacIndexBySymbol(char* symbol) to 
        // return the reaction index for the reaction with char array symbol
        // specifying its symbol.
        // ----------------------------------------------------------------------
        
        static int returnReacIndexBySymbol(char* symbol){
            
            int result;
            for (int i = 0; i < SIZE; i++) {
                result = strcmp(reacLabel[i], symbol);
                if (result == 0){
                    return i;
                }
            }
            return -1;
            
        }
        
        // Static function Utilities::compareTwoCharArrays(char1, char2) 
        // to compare the character arrays char1 and char2. Returns 1
        // (true) if they are equivalent, false (0) if not. Note that
        // if char1 and char2 are of different lengths this will return
        // false. To compare a string to a char array, use 
        // Utilities::stringToChar() to convert the string to a char array.
        // Example usage:
        //
        //    Utilities::compareTwoCharArrays(
        //       Utilities::stringToChar("p+n15-->he4+c12"), reacLabel[21]
        //    );
        //
        // where reacLabel is an array of character arrays.
        
        static bool compareTwoCharArrays(char* char1, char* char2) {
            
            bool result;

            if(strcmp(char1, char2) == 0) {
                result = true;
            } else {
                result = false;
            }
            
//             printf("\nCompare: \"%s\" and %s result=%d", 
//                 char1, char2, result);
            
            return result;

        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::minimumOf(x,y) to return minimum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int minimumOf(int x,int y) {            // integers
            return (x < y) ? x : y;
        }
        
        static double minimumOf(double x,double y) {   // doubles
            return (x < y) ? x : y;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::maximumOf(x,y) to return maximum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int maximumOf(int x,int y) {            // integers
            return (x > y) ? x : y; 
        }
        
        static double maximumOf(double x,double y){    // doubles
            return (x > y) ? x : y; 
        }
        
        
        // ----------------------------------------------------------------------
        // Using the C++ class string instead of char to handle strings 
        // requires #include <iostream>. The C printf command also won't work 
        // correctly with a string type because it isn't typesafe. See the 
        // discussion at
        //
        //    https://stackoverflow.com/questions/10865957/printf-with-stdstring
        //
        // Instead, print with the cout command (which requires #include <iostream>)
        // Example: 
        //
        //      string test("howdy");
        //      string test2(" do");
        //      cout << "Your string is " << test + test2;
        //
        // If instead you try
        //
        //      printf("\n%s",test+test2);
        //
        // it will likely print garbage. However, you can print a string with printf
        // if it is first converted to a char array:
        //
        //      string s = "Howdy World!";
        //      char cs[s.size() + 1];
        //      strcpy(cs,&s[0]);	    // or strcpy(cs,s.c_str());
        //      printf("\n\nstring=%s\n", strcpy(cs,&s[0]));
        //
        // The function Utilities::stringToChar(string) defined below converts a 
        // string to a corresponding character array (returning a pointer to the
        // character array), which either printf or cout can print. Assumes that 
        // the string argument has no more than Utilities::MAXCSIZE characters. Change 
        // Utilities::MAXCSIZE to increase that.
        // ----------------------------------------------------------------------
        
        static char* stringToChar(string s){
            
            // First ensure that string passed to function is not too long
            
            if(s.length() > Utilities::MAXCSIZE - 1){
                
                printf("\n\n*** EXIT: The string\n\n");
                cout << s;
                printf("\nof length %d", s.length());
                printf(" is too long for the Utilities::stringToChar(string) function.");
                printf(
                "\nChange Utilities::MAXCSIZE = %d to value of at least %d and recompile.\n\n", 
                Utilities::MAXCSIZE, s.length() + 1);
                exit(1);
                
            }
            
            static char cs[MAXCSIZE];
            strcpy(cs, &s[0]);          // alternatively strcpy(cs,s.c_str());
            return cs;
        }
        
        // The static function Utilities::charArrayToString (char* a, int size) takes a
        // char array of length size and returns a corresponding string.
        
        static string charArrayToString (char* a, int size) {

            string s = "";
            
            for (int i = 0; i < size; i++) {
                s = s + a[i];
            }
            
            return s;
        }
        
        // ----------------------------------------------------------------------
        // Static function Utilities::startTimer() to start a timer.  Stop timer
        // and display elapsed time with Utilities::stopTimer().
        // ----------------------------------------------------------------------
        
        static void startTimer(){
            
            START_CPU     // Start a timer for rate calculation
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::stopTimer() to stop timer started with 
        // Utilities::startTimer() and display results.
        // ----------------------------------------------------------------------
        
        static void stopTimer(){
            
            STOP_CPU;        // Stop the timer
            printf("\n");
            PRINT_CPU;       // Print timing information for rate calculation
            printf("\n");
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::testTimerCPU() to test CPU timer by 
        // executing long,pointless loop.
        // ----------------------------------------------------------------------
        
        static void testTimerCPU(){
            
            double a,b;
            
            START_CPU;
            for (long count = 1l; count < 500000l; count++) {
                a = sqrt(count);
                b = 1.0/logf(a);
                a = logf(b)/sqrt(a);
            }
            STOP_CPU;
            PRINT_CPU_TEST;
            printf("\n");
        }
    
};    // End class Utilities



/* Class Species to describe all species in the network.  Create new instance
 * for each isotope in the network. Inherits from class Utilities. */

class Species: public Utilities {
    
    // Make data fields private,with external access to them through public setter 
    // and getter functions
    
    private:
        
        int isoindex;        // Index for species in isotope vector
        char IsoLabel[5];    // symbol for isotope
        int ZZ;              // proton number
        int NN;              // neutron number
        int A;               // atomic mass number
        double YY0;          // Abundance Y at beginning of timestep
        double YY;           // current abundance  Y = X/A
        double XX;           // mass fraction  X = Y*A
        double MassExcess;   // mass excess
        double pf[PF];       // partition function entries
        double fplus;        // Total flux currently adding to abundance of this isotope
        double fminus;       // Total flux currently adding to abundance of this isotope
        double keff;         // Effective decay constant = fminus/YY
        double dYdt;         // Current dY/dt for this isotope
        double dXdt;         // Current dX/dt for this isotope
        double pfnow = 1;    // Value of partition function at current T9 for this isotope 
        
    
    public:
        
        // Static method to to renorm all X to sum of one on read-in if
        // invoked after readin of network. Now we start with sumX=1 to machine precision,
        // even if the conversion of the Ys read in to Xs does not give sumX = 1 to machine
        // precision.
        
        static void renormAllX(){
            double oneOver = 1.0/Utilities::sumMassFractions();
            for(int i=0; i<ISOTOPES; i++){
                X[i] *= oneOver;
                Y[i] = X[i]/(double)AA[i]; 
            }
        }
        
        // Public setter functions to set values of private Species class fields
        
        void setisoIndex(int iso){isoindex = iso;}
        
        void setZ(int z){
            ZZ = z;             // Class field
            Z[isoindex] = z;    // Array
        }
        
        void setN(int n){
            NN = n;             // Class field
            N[isoindex] = n;    // Array
        }
        
        void setA(int a){
            A = a;              // Class field
            AA[isoindex] = a;   // Array
        }
        
        void setY(double y){
            YY = y;             // Class field
            XX = (double)A*YY;  // Corresponding X
            Y[isoindex] = y;    // Array
            X[isoindex] = XX;   // Array
        }
        
        void setY0(double y){ YY0 = y; }
        
        void setX(double x){
            XX = x;             // Class field
            YY = XX/(double)A;  // Corresponding Y
            X[isoindex] = x;    // Array
        }
        
        void setM(double m){
            MassExcess = m;             // Class field
            massExcess[isoindex] = m;   // Array
        }
        
        void setisoLabel(char* lpt){
            
            // Fill character array.  *lpt points to initial character.
            // All characters of the character array are accessed by
            // incrementing *lpt.
            
            for(int j=0; j<5; j++){
                IsoLabel[j] = *(lpt+j);               // Class field
                isoLabel[isoindex][j] = IsoLabel[j];  // Array
            }
            
        }
        
        // Set partition function table for this species (24 entries)
        
        void setpf (int i, double pfvalue) { pf[i] = pfvalue; };
        
        // Set PF for this species at current T9
        
        void setpfnow (double value) { pfnow = value; };
        
        void setfminus (double f){
            fminus = f;                 // Class field
            FminusSum[isoindex] = f;    // Array
        }
        
        void setfplus (double f){
            fplus = f;                  // Class field
            FplusSum[isoindex] = f;     // Array entry
        }
        
        void setkeff (double f) { keff=f; }
        
        void setdYdt (double d){
            dYdt = d; 
            dXdt = d*(double)A;
            dYDt[isoindex] = d;
        }
        
        void setdXdt (double d){
            dXdt = d;
            dYdt = d/(double)A;
        }
        
        
        // Public getter functions to return values of class fields
        
        int getZ() {return ZZ; };
        int getN() {return NN; };
        int getA() {return A; };
        double getY0() { return YY0; }
        double getY() {return YY; };
        double getX() {return XX; };
        double getM(){return MassExcess; };
        
        /* Overloaded getLabel function: getLabel() returns pointer to isoLabel character
         array; getLabel(k) returns kth character of isoLabel character array.
         
         Example usage:
        
            Species object1;
              .
              .
              .
            printf("\nLabel=%s",object1.getLabel());
          
            char label[5];
            for(k=0; k<5; k++){
               label[k] = object1.getLabel(k);
            }
            printf("\nLabel=%s",label);
            
         assuming that the data fields of object1 have already been populated using
         the setX functions.
        */
        
        char* getLabel() {return IsoLabel; };         // return pointer to label array
        char getLabel(int k) {return IsoLabel[k]; };  // return kth character of array
        
        double getpf(int i) {return pf[i]; }          // Values in partition function table
        double getpfnow() {return pfnow; }            // PF at current T9 for this isotope
        
        double getfminus(){return fminus; }
        
        double getfplus(){return fplus; }
        
        double getkeff(){ return keff; }
        
        double getdYdt(){return dYdt; }
        
        double getdXdt(){return dXdt; }
        
        
};      // End of class Species



/* Class Reaction to describe all reactions in the network.  Create new instance
 * of this class for each reaction in the network. Inherits from Utilities */

class Reaction: public Utilities {
    
    // Make data fields private,with external access to them through public setter 
    // and getter functions
    
    private:
        
        int reacIndex;               // Index of reaction in current run for each (Z,N)
        int reacClass;               // Reaction class for reaclib (1-8)
        int reacGroupClass;          // Reaction group class (1-6 surrogate for A-F)
        int rgindex;                 // Index of RG containing reaction (0,1,... #RG)
        int RGmemberIndex;           // Index of reaction within its reaction group
        
        char reacGroupClassLett;     // Letter equivalent (A-E) for reacGroupClass
        char* reacGroupSymbol;       // Schematic equil reaction (e.g. a+b<->c)
        int numberReactants;         // Number species on the left side of reaction
        int numberProducts;          // Number species on the right side of reaction
        int numberIsotopes;          // numberReactants + numberProducts in reaction
        char* reacString;            // Character string labeling reaction (not used)
        char resonanceType;          // Resonant (r) or non-resonant (nr) [not used]
        int isEC;                    // Electron capture reaction (1) or not (0)
        int isReverse;               // Reverse reaction (1) or not (0)
        int ispeforward;             // Reactions is labeled "forward" in PE scheme
        bool isEquil;                // In a RG satisfying PE conditions if true
        
        double Q;                    // Q-value for reaction
        double prefac;               // The eta prefac for rates
        double p[7];                 // ReacLib parameters
        int reactantZ[3];            // Array holding Z of reactants
        int reactantN[3];            // Array holding N of reactants
        int reactantA[3];            // A = Z+N
        int productZ[4];             // Array holding Z of products
        int productN[4];             // Array holding N of products
        int productA[4];             // A = Z+N
        int reactantIndex[3];        // Index of species isotope vector for each reactant isotope
        int productIndex[4];         // Index of species isotope vector for each product isotope
        int isoIndex[7];             // Index of species isotope vector for all isotopes in reaction
        
        // Precomputed temperature factors for ReacLib rates.  Computed in 
        // computeConstantFacs(T9,rho),where T9 is the temperature in units of 
        // 10^9 K and rho is the density in g/cm^3.
        
        double T93;                  // T9^{1/3}
        double t1;                   // 1/T9
        double t2;                   // 1/T9^{1/3}
        double t3;                   // T93
        double t4;                   // T9
        double t5;                   // T9^{5/3} = (T93)^5
        double t6;                   // ln (T9)
        
        double Dens[3];              // Powers of density
        double densfac;              // prefac x powers of density
        double rate;                 // Current T-dependent part of rate for reaction
        double Rrate;                // Current full rate (s^-1) for reaction
        double flux;                 // Current effective net flux of reaction if PE
        double flux_true;            // Current true net flux of reaction
        double dErate;               // Current rate of energy release
        char cs[20];                 // Utility character array
        char ccs[20];                // Utility character array
        char* ss;                    // Utility character string
  
  
    public:
        
        // Constructor executed when objects are instantiated
        
        Reaction(){
            
            // Set all reaction objects to not-equilibrated initially.
            
            isEquil = false;
            reacIsActive[reacIndex] =  true;

        }
        
        // Public Reaction setter functions to set values of private class fields
        
        void setreacIndex(int index){reacIndex = index; }
        
        void setreacClass(int index){reacClass = index; }
        
        void setreacGroupClass(int index){
            
            reacGroupClass = index;         // Field in this class
            RGclass[reacIndex] =  index;    // In main
            
            // Use this index to set letters in reacGroupClassLett
            // and symbols in reacGroupSymbol
            
            switch(index){
                
                case 1:
                    reacGroupClassLett = 'A';
                    reacGroupSymbol = Utilities::stringToChar("a<->b");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                    
                case 2:
                    reacGroupClassLett = 'B';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                    
                case 3:
                    reacGroupClassLett = 'C';
                    reacGroupSymbol = Utilities::stringToChar("a+b+c<->d");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                    
                case 4:
                    reacGroupClassLett = 'D';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c+d");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                    
                case 5:
                    reacGroupClassLett = 'E';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c+d+e");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                    
                case 6:
                    reacGroupClassLett = 'F';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c+d+e+f");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
            }
            
        }
        
        void setreacGroupSymbol(char* s){reacGroupSymbol = s; }
        
        void setRGmemberIndex(int i){
            
            RGmemberIndex = i;                   // Field in this class
            RGMemberIndex[reacIndex] = i;        // Write to corresponding array 
            
        }
        
        void setrgindex(int i){rgindex = i; }
        
        void setnumberReactants(int i){
            
            numberReactants = i;                 // Field in this class
            NumReactingSpecies[reacIndex] = i;   // Array
            
        }
        
        void setnumberProducts(int i){
            
            numberProducts = i;                  // Field in this class
            NumProducts[reacIndex] = i;          // Array
            
        }
        
        void setreacString(char* s){ 
            
            // Set field of this class
            
            reacString = s; 
            
            // Set corresponding character array reacLabel 
            
            char p[LABELSIZE];  
            for (int i = 0; i < sizeof(p); i++) { 
                p[i] = s[i]; 
                reacLabel[reacIndex][i] = p[i];
            }
        }
        
        void setisEC(int i){ isEC = i; }
        
        void setisReverse(int i){ isReverse = i; }
        
        void setispeforward(int i){
            ispeforward = i;
            isPEforward[reacIndex] = i;
        }
        
        void setisEquil(bool e){isEquil = e;}
        
        void setQ(double q){ Q = q; }
        
        void setprefac(double p){ prefac = p; }
        
        void setp(double P[]){
            for(int i=0; i<7; i++){
                p[i] = P[i];
            }
        }
        
        void setreactantZ(int z[]){
            for (int i=0; i<numberReactants; i++){
                
                // Change field in this class (Reaction)
                
                reactantZ[i] = z[i];
                
                // Change array 
                
                reacZ[getreacIndex()][i] = z[i];
            }

        }
        
        void setreactantN(int n[]){
            for (int i=0; i<numberReactants; i++){
                
                // Change field in this class (Reaction)
                
                reactantN[i] = n[i];
                
                // Change array
                
                reacN[getreacIndex()][i] = n[i];
            }
        }
        
        void setreactantA(void){
            for (int i=0; i<numberReactants; i++){
                reactantA[i] = reactantZ[i] + reactantN[i];
            }
        }
        
        void setproductZ(int z[]){
            for (int i=0; i<numberProducts; i++){
                
                // Change field in this class (Reaction)
                
                productZ[i] = z[i];
                
                // Change array
                
                prodZ[getreacIndex()][i] = z[i];
            }
        }
        
        void setproductN(int n[]){
            for (int i=0; i<numberProducts; i++){
                
                // Change field in this class (Reaction)
                
                productN[i] = n[i];
                
                // Change array
                
                prodN[getreacIndex()][i] = n[i];
            }
        }
        
        void setproductA(void){
            for (int i=0; i<numberProducts; i++){
                productA[i] = productZ[i] + productN[i];
            }
        }
        
        void setreactantIndex(int n[]){
            for (int i=0; i<numberReactants; i++){
                reactantIndex[i] = n[i];              // Class field
                ReactantIndex[reacIndex][i] = n[i];   // Array
            }
        }
        
        void setproductIndex(int n[]){
            for (int i=0; i<numberProducts; i++){
                productIndex[i] = n[i];               // Class field
                ProductIndex[reacIndex][i] = n[i];    // Array
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes no arguments
        // and constructs isoIndex[] as the concatenation of reactantIndex[]
        // and productIndex[],assuming that those fields have been populated.
        
        void setisoIndex(void){
            for(int i=0; i<numberReactants; i++){
                isoIndex[i] = reactantIndex[i];
            }
            for(int i=0; i<numberProducts; i++){
                isoIndex[i+numberReactants] = productIndex[i];
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes two 
        // arguments and sets a value of a particular isoIndex[].
        
        void setisoIndex(int i,int j){isoIndex[i] =  j;}
        
        void setdensfac(double d){ densfac = d;}
        
        void setrate(double r){ rate = r; }
        
        void setRrate(double r){ Rrate = r; }
        
        void setflux(double f){ flux = f; }
        
        void setflux_true(double f){ flux_true = f;}
        
        // Overloaded dErate setter
        
        void setdErate(double r){ dErate = r; }
        void setdErate(void){ dErate = flux_true*Q; }
        
        
        // Public Reaction getter functions to get values in private fields
        
        int getreacIndex(){ return reacIndex; }

        int getreacClass(){ return reacClass; }
        
        int getreacGroupClass(){ return reacGroupClass; }
        
        // Return reacGroupSymbol as string
        
        char* getreacGroupSymbol(){ return reacGroupSymbol; }
        
        int getRGmemberIndex(){return RGmemberIndex;}
        
        int getrgindex(){return rgindex;}
        
        char getreacGroupClassLett(){ return reacGroupClassLett; }
        
        int getnumberReactants(){ return numberReactants; }
        
        int getnumberProducts(){ return numberProducts; }
        
        // return reacString as string
        
        char* getreacString(){ return reacString; }
        
        // The function getreacChar() returns the string reacString as a
        // pointer to a character array that will work in printf. Alternatively,
        // Utilities::stringToChar() will do same thing.
        
        char* getreacChar(){
            return ccs;
        }
        
        int getisEC(){ return isEC; }
        
        int getisReverse(){ return isReverse; }
        
        int getispeforward(){return ispeforward;}
        
        bool getisEquil(){return isEquil;}
        
        double getQ(){ return Q; }
        
        double getprefac(){ return prefac; }
        
        double getp(int k){ return p[k]; }
        
        int getreactantZ(int k){
            
            if(k > numberReactants-1){
                printf("\n\nERROR: k=%d larger than numberReactants-1 %d; Stopping\n",
                    k,numberReactants);
                exit(-1);                   // Exit program since something is corrupt
            } else {
                return reactantZ[k];
            }
        }
        
        int getreactantN(int k){
            if(k > numberReactants-1){
                printf("\n\nERROR: k-1=%d larger than number reactants %d",
                    k,numberReactants);
                return -1;                 // Exit program since something is corrupt
            } else {
                return reactantN[k];
            }
        }
        
        int getreactantA(int i){
           return (reactantZ[i] + reactantN[i]);
        }
        
        int getproductZ(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k=%d larger than number products %d",
                    k,numberProducts);
                return -1;
            } else {
                return productZ[k];
            }
        }
        
        int getproductN(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d",
                    k,numberProducts);
                return -1;
            } else {
                return productN[k];
            }
        }
        
        int getproductA(int i){
            return (productZ[i] + productN[i]);
        }
        
        int getreactantIndex(int k){
            if(k > numberReactants-1){
                printf("\nReac=%d %s",reacIndex,reacLabel[reacIndex]);
                ss = Utilities::stringToChar(
                    "\nERROR Reaction::getreactantIndex(k): k = %d larger than #reactants-1 = %d");
                printf(stringToChar(ss),k,numberReactants-1);
                return -1;
            } else {
                return reactantIndex[k];
            }
        }
        
        int getproductIndex(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d",
                    k,numberProducts);
                return -1;
            } else {
                return productIndex[k];
            }
        }
        
        double getdensfac(){ return densfac; }
        
        int getisoIndex(int i){return isoIndex[i];}
        
        double getrate(){ return rate; }
        
        double getRrate(){ return Rrate; }
        
        double getflux(){ return flux; }
        
        double getflux_true(){ return flux_true;}
        
        double getdErate(){ return dErate; }
        
        
        // Function Reaction::setupFplusFminus() to set up F+ and F- index for each
        // isotope and to find non-vanishing F+ and F- source terms in network.
        // Function declared static so it can be called as Reaction::setupFplusFminus() 
        // without having to instantiate.
        
        static void setupFplusFminus(){
           
            FplusIsotopeCut[0] = numFluxPlus[0];
            FminusIsotopeCut[0] = numFluxMinus[0];
            for(int i=1; i<numberSpecies; i++){
                FplusIsotopeCut[i] = numFluxPlus[i] + FplusIsotopeCut[i-1];
                FminusIsotopeCut[i] = numFluxMinus[i] + FminusIsotopeCut[i-1];
            }
            
            int currentIso = 0;
            for(int i=0; i<totalFplus; i++){
                FplusIsotopeIndex[i] = currentIso;
                if(i == (FplusIsotopeCut[currentIso]-1)) currentIso ++;
            }
            
            currentIso = 0;
            for(int i=0; i<totalFminus; i++){
                FminusIsotopeIndex[i] = currentIso;
                if(i == (FminusIsotopeCut[currentIso]-1)) currentIso ++;
            }
            
            if(showDetails) fprintf(pFileD,"\n\nMapFplus:\n");
            
            for(int i=0; i<totalFplus; i++){
                MapFplus[i] = tempInt1[i];
                if(showDetails) fprintf(pFileD,"\ni=%d MapFplus[i]=%d  %s",
                    i,MapFplus[i],reacLabel[MapFplus[i]]);
            }
            
            if(showDetails) fprintf(pFileD,"\n\nMapFminus:\n");
            
            for(int i=0; i<totalFminus; i++){
                MapFminus[i] = tempInt2[i];
                if(showDetails) fprintf(pFileD,"\ni=%d MapFminus[i]=%d  %s",
                        i,MapFminus[i],reacLabel[MapFminus[i]]);
            }
            
            // Populate the FplusMin and FplusMax arrays
            
            FplusMin[0] = 0;
            FplusMax[0] = numFluxPlus[0]-1;
            
            for(int i=1; i<numberSpecies; i++){
                FplusMin[i] = FplusMax[i-1] + 1;
                FplusMax[i] = FplusMin[i] + numFluxPlus[i] -1 ;	
            }
            
            // Populate the FminusMin and FminusMax arrays
            
            FminusMin[0] = 0;
            FminusMax[0] = numFluxMinus[0]-1;
            for(int i=1; i<numberSpecies; i++){
                FminusMin[i] = FminusMax[i-1] + 1;
                FminusMax[i] = FminusMin[i] + numFluxMinus[i] -1 ;	
            }
            
            // Populate the FplusFac and FminusFac arrays that hold the factors counting the
            // number of occurrences of the species in the reaction.  Note that this can only
            // be done after ReactionVector::parseF() has been run to give reacMask[i][j].
            
            int tempCountPlus = 0;
            int tempCountMinus = 0;
            for(int i=0; i<ISOTOPES; i++){
                
                for(int j=0; j<SIZE; j++) {
                    if(reacMask[i][j] > 0){
                        FplusFac[tempCountPlus] = (double) reacMask[i][j];
                        tempCountPlus ++;
                    }
                    else if(reacMask[i][j] < 0){
                        FminusFac[tempCountMinus] = (double) -reacMask[i][j];
                        tempCountMinus ++;
                    }	
                }
            }
        }
        
        
        // Static function Reaction::populateFplusFminus() to populate F+ and F- for each
        // isotope set up in setupFplusFminus() from master flux array. Function declared
        // static so it can be called as Reaction::populateFplusFminus() without having
        // to instantiate.
        
        static void populateFplusFminus(){
            
            // Populate the F+ and F- arrays from the master Flux array
            
            for(int i=0; i<totalFplus; i++){
                int indy = MapFplus[i];
                Fplus[i] = FplusFac[i]*Flux[indy];
            }
            
            for(int i=0; i<totalFminus; i++){
                int indy = MapFminus[i];
                Fminus[i] = FminusFac[i]*Flux[indy];
            }
        }
        
        
        // Reaction::computeConstantFacs() to compute the reaction factors 
        // that are constant for a given temperature and density.  Reset each 
        // time T9 or rho changes but stays constant as long as they don't change.  
        // This is required at the beginning of each network integration of the 
        // hydro timestep,since the density and temperature (in addition to
        // abundances because of advection) will generally change over a hydro timestep 
        // in each zone.
        
        void computeConstantFacs(double T9,double rho){

            // Temperature factors in ReacLib rate formula.

            T93 = pow(T9,THIRD); 
            t1 = 1/T9;
            t2 = 1/T93;
            t3 = T93;
            t4 = T9;
            t5 = T93*T93*T93*T93*T93;
            t6 = log(T9);
            
            // Multiply the statistical prefactor by the appropriate 
            // density factors (1 for 1-body,rho for 2-body,and rho^2 
            // for 3-body reactions).
            
            Dens[0] = 1.0f;
            Dens[1] = rho;
            Dens[2] = rho*rho;
            densfac = prefac * Dens[numberReactants - 1];
            setdensfac(densfac); 
            
        }
        
        // Reaction::computeRate(double,double) to compute rates at T and rho. 
        // The quantity rate is the temperature-dependent part,including a possible 
        // partition-function correction.  The quantity Rrate is rate multiplied by
        // appropriate density and statistical factors, which give units of s^-1.  The 
        // flux follows from multiplying Rrate by appropriate abundances Y in computeFlux().
        
        void computeRate(double T9,double rho){

            // The ReacLib library in use has parameters that have been fitted
            // over a range T=1e7 K to T=1e10 K. Thus,it is unreliable to use it
            // for temperatures outside that range. For T < 1e7 K,we will
            // set the corresponding rates equal to zero and continue,with a
            // warning message posted at the end of the calculation that some
            // rates have been set to zero because the temperature was out of bounds.
            // For T > 1e10 K we will warn that the temperature is out of the
            // reliable bounds of the library and exit the program at that point.
            
            if(T9 > 10){           // Temperature above upper bound for library
                
                string str("\n\n\n*** HALTING: t=%7.4es and T9=%5.3f; ");
                str += "rates for T9 > 10 are unreliable for this library.\n\n\n";
                printf(Utilities::stringToChar(str),t,T9);
                exit(1);
                
            }
            
            // If T < 1e7 K set rate = 0 and continue, but set a flag to display
            // a warning message at the end of the calculation. 
            
            if(T9 < 0.01){         // Temperature below lower bound for library
                
                if(!reacLibWarn){
                    
                    reacLibWarn = true;
                    warnTime = t;
                    warnT = T9;
                    
                }
                
                rate = 0.0; 
                
            } else {               // Temperature in bounds for library
                
                rate = exp( 
                      p[0] 
                    + t1*p[1] 
                    + t2*p[2] 
                    + t3*p[3] 
                    + t4*p[4] 
                    + t5*p[5] 
                    + t6*p[6] 
                );
            }
            
            // Correct rate by multiplying by partition function factors for select 
            // reactions if the temperature is high enough.
            
            if(dopf && T9 > pfCut9) pfUpdate();

            // Full rate factor in units of time^-1 (rate from above multiplied 
            // by density factors)
            
            Rrate = densfac*rate;            // Field of this (Reaction) object
            Rate[reacIndex] = Rrate;         // Master rate array
            
        }
        
        // Function Reaction::pfUpdate() to correct the rates using partition 
        // function factors if appropriate.
        
        void pfUpdate(){
            
            double pfnum;
            double pfden;
            double pfFactor;
            
            // Make a partition function correction if this is reverse reaction in
            // sense defined in ReacLib (defined by field Reaction::isReverse = true). 
            // Realistic calculations at higher temperatures should use the partition 
            // function correction so generally set dopf = true unless testing.
            // Partition functions are very near 1.000 if T9 < 1, so we will typically
            // only implement partition function corrections if T9 > pfCut9 = 1.0, but
            // the table of partition functions read in allows pfCut9 as small as 0.1.
            // Because for the temperatures of interest the partition functions for 
            // all light ions (protons, neutrons, alphas, tritons) are equal to 1.0, 
            // the structure of the 8 Reaclib reaction classes specified by 
            // Reaction::reacClass means that this correction is only required for 
            // reverse reactions in Reaclib classes reacClass = 2 and reacClass = 5.
            // Note that the partition function table only extends to T9 = 10, so
            // the code will exit with an error message if you try to calculate 
            // partition functions for T9 > 10.
            
            if(dopf && T9 > pfCut9 && isReverse){
                
                int rin;
                int pin;
                
                if(reacClass == 2){
                    
                    rin = reactantIndex[0];
                    pin = productIndex[1];
                    pfden = currentPF[rin];
                    pfnum = currentPF[pin];
                    
                } else if(reacClass == 5){
                    
                    rin = reactantIndex[1];
                    pin = productIndex[1];
                    pfden = currentPF[rin];
                    pfnum = currentPF[pin];
                    
                } else {
                    
                    pfden = 1.0;
                    pfnum = 1.0;
                }
                
                pfFactor = pfnum/pfden; 
                rate *= pfFactor;

            }
        }
        
        // Function Reaction::showRates() to display computed rates for this
        // Reaction object.
        
        void showRates(){
            if(showDetails) fprintf(pFileD,"\n%d %19s RG=%d densfac=%6.3e rate= %8.5e Rrate=%8.5e",
                getreacIndex(),getreacChar(),getreacGroupClass(),getdensfac(),
                getrate(),getRrate()
            );
        }
        
        
        // Function Reaction::computeFlux() to compute fluxes for reactions.  This is
        // where net flux for a reaction group is set to zero if reaction group is 
        // in equilibrium. The computed flux will be set in the flux field of the
        // Reaction object,and also in the array Flux[].
        
        void computeFlux(){
            
            // If this reaction is in a RG that is not in equilibrium,we
            // need to compute its flux.  The formula for the flux will depend on
            // whether the reaction is 1-body,2-body,or 3-body,and is selected by
            // the switch statement in fluxChooser(numberReactants).
            
            isEquil = !reacIsActive[reacIndex];  // Set isEquil in this Reaction object
            
            if( reacIsActive[reacIndex] ) {
                
                fluxChooser(numberReactants);
            
            // Otherwise the reaction is in an equilibrated reaction group,so set its
            // flux identically to zero.
            
            } else {
                
                fluxChooser(numberReactants);
                flux = 0.0;                 // Put in present (Reaction) object flux field
                Flux[reacIndex] = 0.0;      // Put in main flux array

            }
            
        }   // End of function computeFlux()

  
  
    // Function Reaction::fluxChooser(int) to switch among flux formulas for
    //  one-body,two-body,or three-body reactions.
      
    void fluxChooser(int bodies){

        double kfac;
            
        switch(bodies){
            
            case 1:    // 1-body reactions
                
                kfac = Rrate;
                flux = kfac*Y[ reactantIndex[0] ];	 // In Reaction object flux field
                Flux[reacIndex] = flux;              // In main flux array
                fastSlowRates(kfac);
                flux_true = flux;
            
                break;
                
            case 2:	   // 2-body reactions
                
                kfac = Rrate * Y[ reactantIndex[0] ];
                flux = kfac * Y[ reactantIndex[1] ];  // Put in Reaction object flux field
                Flux[reacIndex] = flux;               // Put in main flux array
                fastSlowRates(kfac);
                flux_true = flux;
                
                break;
                
            case 3:	   // 3-body reactions
                
                kfac = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ];
                flux = kfac * Y[ reactantIndex[2] ];  // Put in Reaction object flux field
                Flux[reacIndex] = flux;               // Put in main flux array
                fastSlowRates(kfac);
                flux_true = flux;
                
                break;
                
        }  // End switch 
        
        //flux_true = flux;

    }       // end function Reaction::fluxChooser(int)
      
      
        // Reaction::fastSlowRates(double) to store fastest and slowest rates 
        // in this timestep,and overall in the calculation. These rates are the
        // rates kfac computed in computeFlux().
        
        void fastSlowRates(double testRate){

            if (isinf(testRate)){
                printf("\n\n***STOP: fastSlowRates() arg infinite for t=%6.4e rindex=%d\n\n",
                    t,reacIndex);
                exit(1);
            }
            
            if (testRate > fastestCurrentRate) {
                fastestCurrentRate = testRate;
                fastestCurrentRateIndex = getreacIndex();
            }
            
            if (testRate < slowestCurrentRate && testRate > 0.0) {
                slowestCurrentRate = testRate;
                slowestCurrentRateIndex = getreacIndex();
            }
            
            if (testRate > fastestOverallRate) {
                fastestOverallRate = testRate;
                fastestOverallRateIndex = getreacIndex();
                timeMaxRate = t;
            }
            
            if (testRate < slowestOverallRate){
                slowestOverallRate = testRate;
                slowestOverallRateIndex = getreacIndex();
            }
            
        }   // End of function Reaction::fastSlowRates()
                
};  // End class Reaction



// Class ReactionVector with  functions that create reaction vectors. Inherits from 
// class Utilities. Uses GSL to define vectors and GSL_BLAS API to manipulate them.

class ReactionVector:  public Utilities {
    
    // Make data fields private,with external access to them through public setter 
    // and getter functions
    
    private:
        
    public:
        
        /* Use the information constructed in parseF (in particular reacMask[j][k]) to create the 
        reaction vectors for the network using the static function 
        ReactionVector::makeReactionVectors(). */
        
        static void makeReactionVectors(){
            
            // Write out the array containing components of the reaction vectors
            
            int uppity = minimumOf(25,numberSpecies);  // limit printout width to 25 species
            if(showDetails) fprintf(pFileD,
                "\n\nREACTION VECTOR ARRAY (%d Reaction vectors with %d species components):\n",
                numberReactions,ISOTOPES);
            if(showDetails) fprintf(pFileD,"\nReaction \\ Species           ");
            for(int k=0; k<uppity; k++){
                if(showDetails) fprintf(pFileD,"%4s  ",isoLabel[k]);
            }
            if(showDetails) fprintf(pFileD,"\n");
            for(int j=0; j<numberReactions; j++){
                if(showDetails) fprintf(pFileD,"%4d %22s [ ",j,reacLabel[j]);
                for(int k=0; k<uppity-1; k++){
                    if(showDetails) fprintf(pFileD,"%2d    ", reacMask[k][j]);
                }
                if(showDetails) fprintf(pFileD,"%2d ]\n", reacMask[uppity-1][j]);
            }
            
            
            // Fill std::vector component entries with data contained in  
            // reacMask[j][i] (notice reversed indices because outer loop is reactions
            // and inner loop is isotopes and reacMask is defined as reacMask[iso][reac]
            // where iso in the isotopic index and reac is the reac index.)
            
            for (int i = 0; i < SIZE; i++) {     // Loop over std::vector array RV

                // Populate RV[i] with components
                
                for(int j=0; j < ISOTOPES; j++){
                    
                    RV[i].push_back(reacMask[j][i]);

                }

            }
            
            // Display std:vector components just created
            
            if(showDetails2) fprintf(pfnet,"\n\nSTD REACTION VECTORS\n");
            
            // Outer loop over elements of vector array
            
            for (int i = 0; i < SIZE; i++) {

                if(showDetails2) fprintf(pfnet, "\nRV[%d]: [", i);
                
                // Inner iterates over components of each vector. Starting iterator 
                // is specified by .begin() and ending iterator is specified by 
                // .end(). The auto keyword asks the compiler to deduce the type 
                // of the variable iter from the context.
                
                for (auto iter = RV[i].begin();
                     
                     iter != RV[i].end(); iter++) {
                    
                    // *iter gets the value pointed to by the iterator
                    
                    if(showDetails2) fprintf(pfnet, "%3d", *iter);
                
                }
                    if(showDetails2) fprintf(pfnet, " ]");
            }
            
            // Now compare reaction vectors pairwise as check
            
            std::vector<int> vec1;
            std::vector<int> vec2;
            
            if(showDetails2){
            
                for (int i=0; i<SIZE; i++){
                    
                    vec1 =  RV[i];
                
                    for (int j=i; j<SIZE; j++){
                        
                        vec2 =  RV[j];
                        int checker = compareVectors(vec1, vec2);
                        
                        if(checker > 0 && showDetails2)
                            fprintf(pfnet, "\ni=%d j=%d %s <-> %s ck=%d", 
                            i, j, reacLabel[i], reacLabel[j], checker);
                    }
                    
                }
            }
            
        }  // End function makeReactionVectors()
        
        
        /**
         * Compare two std vectors
         * @param rv1 First vector to compare
         * @param rv2 Second vector to compare
         * @return 0 if the vectors are not equal, 1 if they are the same, 2 if one
         * vector is the negative of the other
         */
        
        static int compareVectors(const std::vector<int> & rv1, const std::vector<int> & rv2) {
            
            int retVal = 0;
            
            // Check for element-wise equality
            
            auto equal = std::equal(rv1.begin(),rv1.end(),rv2.begin(),rv2.end());
            
            // If not equal, are they negated images?
            
            if (!equal) {
                
                // Flip rv2 by using the handy transformation instead of a for loop
                
                std::vector<int> flippedRv2(rv2);
                std::transform(rv2.begin(), rv2.end(), flippedRv2.begin(),
                    std::bind(std::multiplies<int>(), std::placeholders::_1, -1.0));
                
                // Check the value and set the flag
                
                bool flipped = std::equal(rv1.begin(), rv1.end(), flippedRv2.begin(), flippedRv2.end());
                if (flipped) retVal = 2;
                
            } else {
                
                retVal = 1;  
            }
            
            return retVal;
        }
        

    
    // ------------------------------------------------------------------------
    // ReactionVector::sortReactionGroups() uses compareVectors to sort all 
    // reactions in the network into reaction groups labeled by a series of 
    // integers 0,1,...  All reactions in a reaction group have the same 
    // reaction vector up to a sign. The array RGindex[] of dimension SIZE 
    // holds the integer labeling reaction group (0,1,... #RG) for each 
    // reaction after this function is executed. Make static since function
    // is called only once in a calculation.
    // ------------------------------------------------------------------------
    
    static void sortReactionGroups(void){
        
        // Cycle over all reaction vectors and compare them pairwise to 
        // assign to reaction groups.  The integer rg labels the reaction 
        // group currently being filled.  The integer ck indicates whether 
        // a pair of vectors are equivalent (ck = 1), are the negative of 
        // each other (ck = 2), or are not equivalent (ck = 0). Requires
        // approximately SIZE*SIZE/2 vector comparisons.
        
        int rg = -1;
        int ck = -1;
        
        // Initialize reaction vector indices that give the reaction group
        // for a reaction vector to -1, so that we can tell in sorting
        // into reaction groups whether a reaction vector has already been 
        // processed.
        
        for(int i=0; i<SIZE; i++){
            RGindex[i] = -1;
        }
        
        int numberMembers;
        
        std::vector<int> vec1;
        std::vector<int> vec2;
        
        int compares = 0;
        
        for (int i=0; i<SIZE; i++){
            
            if(numberMembers > 0) rg ++;
            numberMembers = 0;
            vec1 = RV[i];

            // Since we only need to compare pairwise, once for each pair,
            // inner loop can start at j=i.
            
            for(int j=i; j<SIZE; j++){
                
                compares ++;
                
                vec2 = RV[j]; 
                
                // Compare the vectors vec1 and vec2
                
                ck = compareVectors(vec1, vec2);
                
                int preRGindex = RGindex[j];
                
                if(ck > 0 && RGindex[j] < 0) {
                    
                    RGindex[j] = rg;
                    numberMembers ++;
                }
            }
            
            // Store the number of member reactions in this reaction group 
            // for later use
            
            RGnumberMembers[rg] = numberMembers;

        }
        
        printf("\nVector pair comparisons = %d", compares);
        cout.flush();
        
        // If the last trial reaction group has no members, subtract 
        // one from rg (which was incremented at the beginning of the trial).
        
        if(numberMembers == 0) rg --;
        
        numberRG = rg + 1;   // Total # reaction groups (add 1 since rg starts at 0)

        // Diagnostic showing reaction group associated with each reaction
        
        if(showDetails2) fprintf(pfnet, "\n");
        if(showDetails2) fprintf(pfnet, "\nCHECK: Each reaction should be in one and only one RG\n");
        for(int i=0; i<SIZE; i++){
            if(showDetails2) fprintf(pfnet, "\n%d %s RG=%d", i, reacLabel[i], RGindex[i]);
        }
        
        // Output the components of the reaction groups pfnet -> network.out.
        
        numberSingletRG = 0;
        
        fprintf(pfnet,"\n\n\nPARTIAL EQUILIBRIUM REACTION GROUPS");
        for(int i=0; i<numberRG; i++){
            
            fprintf(pfnet,"\n\nReaction Group %d:", i);
            int rgindex = -1;
            int rcounter = 0;
            
            for(int j=0; j<SIZE; j++){
                if(RGindex[j] == i){
                    rgindex ++; 
                    rcounter ++;
                    setRG(j, RGclass[j], RGindex[j]);
                    fprintf(pfnet,
                    "\n%s reacIndex=%d RGindex=%d RGclass=%d RGreacIndex=%d isForward=%d RG:%s",
                    reacLabel[j],j, rgindex, RGclass[j], RGMemberIndex[j],
                    isPEforward[j], Utilities::stringToChar(RGstring[j]));
                }
            }
            
            // Number of reactions in each RG
            
            RGnumberMembers[i] = rcounter;
            
            // Keep track of RG that are singlets (one reaction in RG)
            
            if (rcounter == 1) numberSingletRG ++;
        }
        
        if(showDetails2) fprintf(pfnet,"\n");
        
        if(showDetails2) fprintf(pfnet, "\n");
        if(showDetails2) fprintf(pfnet, "\nCHECK: Sum over members of ");
        if(showDetails2) fprintf(pfnet, "each RG should equal total number of reactions (SIZE).\n");
        
        int resum = 0;
        for(int i=0; i<numberRG; i++){
            resum = resum + RGnumberMembers[i];
            if(showDetails2) fprintf(pfnet, "\nRG=%d #reactions=%d", i, RGnumberMembers[i]);
        }
        
        if(showDetails2) fprintf(pfnet,"\n\nSum=%d SIZE=%d\n", resum, SIZE);
        
        fflush(pfnet);  // Dump buffer to force quicker print
        
    }   // End function sortReactionGroups()
        
    
        /*
        * ReactionVector::parseF() to find contributions to F+ and F- of each reaction for each 
        * isotope. This is executed only once at the beginning of the calculation to determine 
        * the structure of the network. Since executed only once, make it static so it can
        * be called directly from the class using ReactionVector::parseF(), without having
        * to instantiate.
        */
        
        static void parseF(){
            
            int incrementPlus = 0;
            int incrementMinus = 0;
            
            // Loop over all isotopes in the network
            
            for(int i=0; i<numberSpecies; i++) {
                int total = 0;
                int numFplus = 0;
                int numFminus = 0;
                
                // Loop over all possible reactions for this isotope,finding those that
                // change its population up (contributing to F+) or down (contributing to F-).
                
                for(int j=0; j<numberReactions; j++) {
                    
                    int totalL = 0;
                    int totalR = 0;
                    
                    // Loop over reactants for this reaction
                    
                    for(int k=0; k<NumReactingSpecies[j]; k++) {
                        if(Z[i] == reacZ[j][k] && N[i] == reacN[j][k]) totalL ++;
                    }
                    
                    // Loop over products for this reaction
                    
                    for(int k=0; k<NumProducts[j]; k++) {
                        if(Z[i] == prodZ[j][k] && N[i] == prodN[j][k]) totalR ++;
                    }
                    
                    total = totalL - totalR;
                    
                    if(total > 0){       // Contributes to F- for this isotope
                        numFminus ++;
                        reacMask[i][j] =  -total;
                        tempInt2[incrementMinus + numFminus-1] = j;
                    } 
                    else if(total < 0){  // Contributes to F+ for this isotope
                        numFplus ++;
                        reacMask[i][j] =  -total;
                        tempInt1[incrementPlus + numFplus-1] = j;
                    } else {             // Doesn't contribute to flux for this isotope
                        reacMask[i][j] = 0;
                    }
                }
                
                // Keep track of the total number of F+ and F- terms in the network 
                // for all isotopes
                
                totalFplus += numFplus;
                totalFminus += numFminus;
                
                numFluxPlus[i] = numFplus;
                numFluxMinus[i] = numFminus;
                
                incrementPlus += numFplus;
                incrementMinus += numFminus;
                
            }
            
            // Display isotope component array
            
            if(showDetails) fprintf(pFileD,
            "\n\n\nFLUX-ISOTOPE COMPONENT ARRAY (negative n for F-; positive n for F+ for given isotope):");
            if(showDetails) fprintf(pFileD,"\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
            
            int uppity = minimumOf(30,numberSpecies);  // limit printout width to 30 species
            if(showDetails) fprintf(pFileD,"\n\nIndex             Reaction");
            for(int k=0; k<uppity; k++){
                if(showDetails) fprintf(pFileD,"%5s",isoLabel[k]);
            }
            for(int j=0; j<numberReactions; j++){
                if(showDetails) fprintf(pFileD,"\n%3d %22s",j,reacLabel[j]);
                for(int k=0; k<uppity; k++){
                    if(showDetails) fprintf(pFileD," %4d", reacMask[k][j]);
                }
            }
            
            if(showDetails) fprintf(pFileD,
                "\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d,out of %d x %d = %d possibilities.\n",
                totalFplus,totalFminus,SIZE,ISOTOPES,SIZE*ISOTOPES);
            
        }   // End of function parseF()
    
    
};    // end class ReactionVector



/*
*   class MatrixUtils inherits from Utilities. Functions to create GSL matrices and vectors
*   given standard C++ arrays and to compute matrix-vector multiply M*v using BLAS.
*/


class MatrixUtils: public Utilities {

    private:

        // Rows and columns in flux matrix.  Number of columns in flux matrix is 
        // always equal to rows in abundance vector.
        
        const static int FLUXROWS = 300;
        const static int FLUXCOLS = 100;

    public:

        // Allocate and populate GSL abundance vector
        
        void buildGSLVector(double a[]){
            abundances = gsl_vector_alloc(FLUXCOLS);

            for (int i = 0; i < FLUXCOLS; i++){
                gsl_vector_set(abundances,i,a[i]);
            }
        }

        // Allocate and populate GSL flux matrix
        
        void buildGSLMatrix(double f [][FLUXCOLS]){
            
            fluxes = gsl_matrix_alloc(FLUXROWS, FLUXCOLS);

            for (int i = 0; i < FLUXCOLS; i++){
                for (int j = 0; j < FLUXROWS; j++){
                    gsl_matrix_set(fluxes, i, j, f[i][j]);
                }
            }
        }

        // Matrix vector multiply fluxes * abundances and store back into abundances
        
        gsl_vector multiply(gsl_vector a, gsl_matrix f){
            
            gsl_blas_dgemv(CblasNoTrans, 1.0, fluxes, abundances, 0.0, abundances);

            return *abundances;
        }

};  // end of class MatrixUtils



// Class ReactionGroup to handle partial equilibrium (PE) reaction groups (RG).  
// Inherits from class Utilities. In this code we instantiate an array of ReactionGroup 
// objects RG[], one for each reaction group in the network. For consistency we assume
// that all reactions in the network are in one (and only one) reaction group. Thus, 
// there will be reaction groups with a single member reaction if the network does not
// have inverse reactions for every reaction.

class ReactionGroup:  public Utilities {
    
    // Make data fields private,with external access to them through public setter 
    // and getter functions.  Static functions can be called directly from the class
    // without having to instantiate.
    
    public:
        
        static const int maxreac = 10;         // Max possible reactions in this RG instance
        int nspecies[6] = {2,3,4,4,5,6};       // Number isotopic species in 6 RG classes
        int niso;                              // Number isotopic species in this RG object
        int RGn;                               // Index this object in RG array (0,1,... #RG)
        int numberMemberReactions;             // Number of reactions in this RG instance
        int memberReactions[maxreac];          // reacIndex of reactions in reaction group
        int numberReactants[maxreac];          // # reactants for each reaction in RG
        int numberProducts[maxreac];           // #products for each reaction in RG
        int refreac = -1;                      // Ref. reaction for this RG in memberReactions

        int rgclass = -1;                      // Reaction group class (1-5)
        bool isEquil;                          // True if RG in equilibrium; false otherwise
        bool isForward[maxreac];               // Whether reaction in RG labeled forward
        double flux[maxreac];                  // Current flux for each reaction in RG
        double netflux;                        // Net flux for the entire reaction group
        char reaclabel[maxreac][LABELSIZE];    // Member reaction label
        
        // Partial equilibrium quantities
        
        double crg[4];                 // Constants c1,c2,... (allocate dynamically?)
        int numberC;                   // Number constants crg[] for this RG object
        double rgkf;                   // Forward rate parameter for partial equilibrium
        double rgkr;                   // Reverse rate parameter for partial equilibrium
        
        double aa,bb,cc;               // Quadratic coefficients a,b,c
        double alpha,beta,gamma;       // Coefficients for cubic ~ quadratic approximation
        double qq;                     // q = 4ac-b^2
        double rootq;                  // Math.sqrt(-q)
        double tau;                    // Timescale for equilibrium

        double equilRatio = 0;         // Equilibrium ratio of abundances
        double kratio = 0;             // Ratio k_r/k_f. Equal to equilRatio at equilibrium
        double eqcheck[5];             // Population ratio to check equilibrium
        double Yminner;                // Current minimum Y in reaction group
        double eqRatio[5];             // Ratio eqcheck[i]/equiTol. 
        double maxRatio;               // Largest eqcheck[]/equiTol in reaction group
        double minRatio;               // Smallest eqcheck[]/equiTol in reaction group 
        double thisDeviousRG;          // Value of thisDevious for this RG

        double lambda;                 // Progress variable for RG (computed; not used)
        double lambdaEq;               // Equil value progress variable (computed; not used)
        
        // Reactions of the Reaclib library have 1-3 reactants and 1-4 products,so create arrays
        // to accomodate all possibilities without having to resize arrays.
        
        int reactantIsoIndex[3];       // Species index of reactants
        int productIsoIndex[4];        // Species index of products
        
        // Reactions in the ReacLib library involve 2-5 species.  Make the following
        // arrays of dimension 5 to accomodate all possibilities without having to
        // resize arrays.
        
        int isoindex[5];               // Species index for participants in reaction   
        char isolabel[5][5] = 
           {' ', ' ', ' ', ' ', ' '};  // Label of isotopic species in RG reactions
        int isoZ[5] = {0,0,0,0,0};     // Z for niso isotopes in the reactions of group
        int isoN[5];                   // N for niso isotopes in the reactions of the group
        double isoA[5];                // A for niso isotopes in the reactions of the group
        double isoYeq[5]={0,0,0,0,0};  // Y_eq for niso isotopes in reactions of the group
        double isoY[5]={0,0,0,0,0};    // Current Y for niso isotopes in reactions of group
        double isoY0[5];               // Y0 for niso isotopes in the reactions of the group

    
    public:
    
    // Constructor executed when ReactionGroup objects instantiated
        
    ReactionGroup(int rgn){
        RGn = rgn;
        isEquil = false;
    }
    
    // Public ReactionGroup setter functions to set values of private class fields
    
    void setnumberMemberReactions(int n){numberMemberReactions = n;}
    
    void setmemberReactions (int i,int index){memberReactions[i] = index;}
    
    void setnumberReactants(int i,int j){numberReactants[i] = j;}
    
    void setnumberProducts(int i,int j){numberProducts[i] = j;}
    
    void setisEquil(bool b) {isEquil = b;}
    
    void setisForward(int i,bool b){isForward[i] = b;}
    
    int setrefreac(){
        
        // Set the reference reaction for the ReactionGroup by looping
        // through and choosing the first reaction that has .forward = true.
        
        for (int i = 0; i < numberMemberReactions; i++) {
            
            if (isPEforward[i]) {
                refreac = i;
                return refreac;
            }
        }
        
        // Take care of anomalous case where there are no forward reactions with 
        // a given reaction vector (which could be generated by suppressing
        // all the forward reactions in a reaction group,for example).
        
        if (refreac == -1) {
            refreac = 0;
            if(showDetails) fprintf(pFileD,"\n*** Reaction group %d has no forward reactions ***",
                RGn);
        }
        
        return refreac;
    }
    
    void setflux(int i,double f){flux[i] = f;}

    void setRGn(int rgn){RGn = rgn;}
    
    void setrgclass(int rc){rgclass = rc;}

    void setreaclabel(int k,char *s){ 
        
        // Convert from string to char array
        
        char p[LABELSIZE];  
        for (int i = 0; i < sizeof(p); i++) { 
            p[i] = s[i]; 
            reaclabel[k][i] = p[i];
        }
    }
       
    // ReactionGroup::setRGfluxes() to set all fluxes in RG from the array Flux[],
    // and compute the net flux for the RG.
    
    void setRGfluxes(){

        for(int i = 0; i < numberMemberReactions; i++){
            setflux(i,Flux[ memberReactions[i] ]);
        }
        
        // Set current net flux for the reaction group object
        
        netflux = sumRGfluxes();
    }
    
    void setnetflux(double f){netflux = f;}
    
    void setniso(int rgclass){niso = nspecies[rgclass-1];}
    
    void setisoindex(int i,int j){isoindex[i] = j;}
    
    void setisolabel(int k,char s[]){ 
        for (int i = 0; i < 5; i++) { 
            isolabel[k][i] = s[i];
        }
    }
    
    void setisoZ(int i,int j){isoZ[i] = j;}
    
    void setisoN(int i,int j){isoN[i] = j;}
    
    void setisoA(int i,int j){isoA[i] = j;}
    
    void setisoY0(int i,double y){isoY0[i] = y;}
    
    void setisoY(int i,double y){isoY[i] = y;}
    
    void setisoYeq(int i,double y){isoYeq[i] = y;}
    
    void setreactantIsoIndex(int i,int j){
        reactantIsoIndex[i] = j;
    }

    void setproductIsoIndex(int i, int j){
        productIsoIndex[i] = j;
    }
    
    void seteqcheck(int k, double eq){eqcheck[k] = eq;}
    void seteqRatio(int k, double eq){eqRatio[k] = eq;}
    
    void setthisDeviousRG(double f){thisDeviousRG = f;}
    
    void setnumberC(int k){numberC = k;}
    
    
    // Public ReactionGroup getter functions to retrieve values of private class fields
    
    int getnumberMemberReactions(){return numberMemberReactions;}
    
    int getmemberReactions (int i) {return memberReactions[i];}
    
    int getnumberReactants(int i){return numberReactants[i];}
    
    int getnumberProducts(int i){return numberProducts[i];}
    
    bool getisEquil() {return isEquil;}
    
    bool getisForward(int i){return isForward[i];}
    
    int getrefreac(){return refreac;}
    
    double getflux(int i){return flux[i];}
    
    int getRGn(){return RGn;}
    
    int getrgclass(){return rgclass;}
    
    char* getreaclabel(int k){return reaclabel[k];}
    
    double getRGfluxes(int index){return flux[index];}
    
    double getnetflux(){return netflux;}
    
    int getniso(){return niso;};
    
    int getisoindex(int i){return isoindex[i];}
    
    char* getisolabel(int k){return isolabel[k];}
    
    int getisoZ(int i){return isoZ[i];}
    
    int getisoN(int i){return isoN[i];}
    
    int getisoA(int i){return isoA[i];}
    
    double getisoY0(int i){return isoY0[i];}
    
    double getisoY(int i){return isoY[i];}
    
    double getisoYeq(int i){return isoYeq[i];}
    
    int getreactantIsoIndex(int i){return reactantIsoIndex[i];}
    
    int getproductIsoIndex(int i){return productIsoIndex[i];}
    
    double geteqcheck(int k){return eqcheck[k];}
    
    double geteqRatio(int k){return eqRatio[k];}
    
    double getequilRatio(){return equilRatio;}
    
    double getthisDeviousRG(){return thisDeviousRG;}
    
    double getmaxRatio(){return maxRatio;}
    
    double getrgkf(){return rgkf;}
    
    double getrgkr(){return rgkr;}
    
    double getkratio(){return kratio;}
    
    double getaa(){return aa;}
    
    double getbb(){return bb;}
    
    double getcc(){return cc;}
    
    double getalpha(){return alpha;}
    
    double getbeta(){return beta;}
    
    double getgamma(){return gamma;}
    
    double getqq(){return qq;}
    
    double gettau(){return tau;}
    
    bool getisEqual(){return isEquil;}
    
    int getnumberC(){return numberC;}
    
    double getcrg(int k){return crg[k];}
    
    
    // Function ReactionGroup::showRGfluxes to show all current fluxes in 
    // reaction groups
    
   void showRGfluxes(){
        
        if(showDetails) fprintf(pFileD, "\n\nRG=%d", RGn);
        
        double fac;
        for(int i=0; i<numberMemberReactions; i++){
            if(isForward[i]){
                fac = 1.0;
            } else {
                fac = -1.0;
            }
            if(showDetails) fprintf(pFileD,
                "\nshowRGfluxes: %d %s RGclass=%d isForward=%d t=%7.4e dt=%7.4e flux=%7.4e",
                i,reacLabel[memberReactions[i]],rgclass,isForward[i],t,dt,
                fac*flux[i] );
            
        }
        
        if(showDetails) fprintf(pFileD,"\n");
        if(isEquil){
            if(showDetails) fprintf(pFileD,"showRGfluxes: NetRGflux=%7.4e\nEQUILIBRATED", netflux); 
        } else {
            if(showDetails) fprintf(pFileD,"showRGfluxes: NetRGflux=%7.4e\nNOT EQUILIBRATED",netflux); 
        }
        
    }
    
    // Function ReactionGroup::sumRGfluxes to sum net flux for this reaction group.
    // This corresponds to sumFluxes in original Java program.
    
    double sumRGfluxes(){
        
        double sumf = 0.0;
        double fac;

        for (int i=0; i<numberMemberReactions; i++){
            fac = -1.0;
            if(isForward[i]) fac = 1.0;
            sumf += fac*flux[i];
        }

        netflux = sumf;
        return sumf;
        
    }    // End function sumRGfluxes()
    
    
    // ---------------------------------------------------------
    // Function ReactionGroup::computeEquilibrium() to compute
    // all partial equilibrium quantities
    //----------------------------------------------------------
    
    void computeEquilibrium() {
        
        mostDevious = 0.0;
        computeEquilibriumRates();
        putY0();
        computeC();
        computeQuad();
        
        // Compute net flux and progress variable for this reaction group
        
        if (isEquil) {
            netflux = sumRGfluxes();
            lambda = netflux * dt;   // lambda not presently used
        }
    }
    
    
    // -----------------------------------------------------------------
    // Function ReactionGroup::computeEquilibriumRates() to compute the 
    // net forward and reverse rates k_f and k_r required in partial 
    // equilibrium approximation.
    // -----------------------------------------------------------------
    
    void computeEquilibriumRates() {
        
        double kf = 0.0;
        double kr = 0.0;
        
        // Sum all contributions from members of reaction group
        
        for (int j = 0; j < numberMemberReactions; j++) {
            int rin = memberReactions[j];
            if (isForward[j]) {
                kf += Rate[rin];
            } else {
                kr += Rate[rin];
            }
        }
        
        // Store forward and reverse rates
        
        rgkf = kf;
        rgkr = kr;
        
    }  // End of function computeEquilibriumRates()
    
    
    // -----------------------------------------------------------------------
    // Function ReactionGroup::putY0() to put the values of Y0 and Y at 
    // beginning of timestep into the Y0[] and Y[] arrays for this 
    // RG object.
    // -----------------------------------------------------------------------
    
    void putY0() {
        
        int ii;
        
        for (int k = 0; k < niso; k++) {
            ii = isoindex[k];
            isoY0[k] = Y0[ii];
            isoY[k] = Y[ii];   // old error: isoY[k]=isoY0[k];
        }
        
    }
    
   
   // -----------------------------------------------------------------------
   // Function ReactionGroup::computeC() to compute values of constants crg[]
   // -----------------------------------------------------------------------
   
    void computeC() {
        
        switch (rgclass) {
            
            // Reaclib class 7,which can't equilibrate for standard 
            // ReacLib library because there are no inverse reactions
            // in the library for the reaclib class.
            
            case -1:
                
                break;
                
            case 1:    // a <-> b
                
                crg[0] = isoY0[0] + isoY0[1];
                
                break;
                
            case 2:    // a+b <-> c
                
                crg[0] = isoY0[1] - isoY0[0];
                crg[1] = isoY0[1] + isoY0[2];
                
                break;
                
            case 3:    // a+b+c <-> d
                
                crg[0] = isoY0[0] - isoY0[1];
                crg[1] = isoY0[0] - isoY0[2];
                crg[2] = THIRD * (isoY0[0] + isoY0[1] + isoY0[2]) + isoY0[3];
                
                break;
                
            case 4:    // a+b <-> c+d
                
                crg[0] = isoY0[0] - isoY0[1];
                crg[1] = isoY0[0] + isoY0[2];
                crg[2] = isoY0[0] + isoY0[3];
                
                break;
                
            case 5:    //  a+b <-> c+d+e
                
                crg[0] = isoY0[0] + THIRD * (isoY0[2] + isoY0[3] + isoY0[4]);
                crg[1] = isoY0[0] - isoY0[1];
                crg[2] = isoY0[2] - isoY0[3];
                crg[3] = isoY0[2] - isoY0[4];
                
                break;
        }
        
    }     // End of function computeC()
    
    
    // --------------------------------------------------------------
    // Function ReactionGroup::computeQuad() to compute the quadratic 
    // coefficients needed for the equilibrium solution and to 
    // compute the equilibrium solution
    // --------------------------------------------------------------
    
    void computeQuad() {
        
        switch (rgclass) {
            
            // Reaclib class 7,which can't equilibrate in the standard
            // ReacLib library because there is no inverse reaction
            // in the reaction library
            
            case -1: 
                
                break;
                
            case 1:  // a <-> b
                
                aa = 0;
                bb = -rgkf;
                cc = rgkr;
                
                break;
                
            case 2:  // a+b <-> c
                
                aa = -rgkf;
                bb = -(crg[0] * rgkf + rgkr);
                cc = rgkr * (crg[1] - crg[0]);
                
                break;
                
            case 3:  // a+b+c <-> d
                
                aa = -rgkf * isoY0[0] + rgkf * (crg[0] + crg[1]);
                bb = -(rgkf * crg[0] * crg[1] + rgkr);
                cc = rgkr * (crg[2] + THIRD * (crg[0] + crg[1]));
                
                break;
                
            case 4:  // a+b <-> c+d
                
                aa = rgkr - rgkf;
                bb = -rgkr * (crg[1] + crg[2]) + rgkf * crg[0];
                cc = rgkr * crg[1] * crg[2];
                
                break;
                
            case 5:  //  a+b <-> c+d+e
                
                alpha = crg[0] + THIRD * (crg[2] + crg[3]);
                beta = crg[0] - TWOTHIRD * crg[2] + THIRD * crg[3];
                gamma = crg[0] + THIRD * crg[2] - TWOTHIRD * crg[3];
                aa = (3 * crg[0] - isoY0[0]) * rgkr - rgkf;
                bb = crg[1] * rgkf
                - (alpha * beta + alpha * gamma + beta * gamma) * rgkr;
                cc = rgkr * alpha * beta * gamma;
                
                break;
        }
        
        // Compute the q = 4ac - b^2 parameter,equil timescale tau,and
        // isoYeq[0] (which is then be used to compute the other isoYeq[].
        
        if (rgclass > 1) {
            qq = computeq(aa,bb,cc);
            rootq = sqrt( max(-qq, 0.0) );
            if (numberMemberReactions > 1) {
                tau = 1.0 / rootq;
            }
            isoYeq[0] = computeYeq(aa,bb,rootq);
            
        } else {
            qq = -1.0;
            tau = 1.0 / rgkf;
            isoYeq[0] = rgkr / rgkf;
        }
        
        // Compute the other equilibrium populations in the reaction pair
        // and abundance ratios
        
        switch (rgclass) {
            
            // Reaclib class 7, which can't equilibrate because no inverse in 
            // standard Reaclib library.
            
            case -1: 
            
                break;
                
            case 1:    // a <-> b
                
                isoYeq[1] = crg[0] - isoYeq[0];
                equilRatio = isoY[0] / isoY[1];
                
                break;
                
            case 2:    // a+b <-> c
                
                isoYeq[1] = crg[0] + isoYeq[0];
                isoYeq[2] = crg[1] - isoYeq[1];
                equilRatio = isoY[0] * isoY[1] / isoY[2];

                break;
                
            case 3:    // a+b+c <-> d
                
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = isoYeq[0] - crg[1];
                isoYeq[3] = crg[2] - isoYeq[0] + THIRD * (crg[0] + crg[1]);
                equilRatio = isoY[0] * isoY[1] * isoY[2] / isoY[3];
                
                break;
                
            case 4:    // a+b <-> c+d
                
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = crg[1] - isoYeq[0];
                isoYeq[3] = crg[2] - isoYeq[0];
                equilRatio = isoY[0] * isoY[1] / ( isoY[2] * isoY[3] );
                
                break;
                
            case 5:    //  a+b <-> c+d+e
                
                isoYeq[1] = isoYeq[0] - crg[1];
                isoYeq[2] = alpha - isoYeq[0];
                isoYeq[3] = beta - isoYeq[0];
                isoYeq[4] = gamma - isoYeq[0];
                equilRatio = isoY[0] * isoY[1] / ( isoY[2] * isoY[3] * isoY[4] );
                
                break;
        }
        
        
        // Compute the equilibrium value of the progress variable
        
        lambdaEq = isoY0[0] - isoYeq[0];         // Not presently used
        
        // Compute the population ratios used to check equilibration
        
        kratio = rgkr / rgkf;
        
        if(t > equilTime) computeEqRatios();
        
    }    // End function computeQuad()
    
    
    // ---------------------------------------------------------------------
    // Function ReactionGroup::computeq to compute q = 4ac-b^2 for quadratic 
    // solution
    // ---------------------------------------------------------------------
    
    double computeq(double a,double b,double c) {
        return 4 * a * c - b * b;
    }
    
    // ---------------------------------------------------------------------
    //  Function ReactionGroup::computeYeq to compute Yeq[0]
    // ---------------------------------------------------------------------
    
    double computeYeq(double a,double b,double rootq) {
        return -0.5 * (b + rootq) / a;
    }
    
    
    // ---------------------------------------------------------------------
    // Function ReactionGroup::computeEqRatios to compute array of population 
    // ratios and check for equilibration in each reaction group.
    // ---------------------------------------------------------------------
    
    void computeEqRatios() {
        
        // Compute thisDevious,which measures difference from equil.
        // Limit how small denominator can be in following to prevent 
        // possible divide by zero.
        
        thisDevious = abs((equilRatio - kratio) / max(kratio, GZ));
        
        thisDeviousRG = thisDevious;  // Store also in RG group field
        
        // Store max value of thisDevious
        
        if (isEquil && thisDevious > mostDevious) {
            mostDevious = thisDevious;                
            mostDeviousIndex = RGn;
        }
        
        // The return statements in the following if-clauses cause reaction
        // groups already in equilibrium to stay in equilibrium. Otherwise, if
        // the RG is in equilibrium (isEquil=true) but the tolerance condition
        // thisDevious < deviousMax is no longer satisfied, the RG is removed 
        // from equilibrium.
        
        if (isEquil && thisDevious < deviousMax) {
            return;
        } else if (isEquil && thisDevious >= deviousMax) {
            removeFromEquilibrium(1);
            return;
        }
        
        Yminner = 1000;
        maxRatio = 0;
        minRatio= 1000;
        
        
        // If we have gotten this far without returning from this function, the
        // reaction group was not in equilibrium before.  See if it is now.
        // Determine if reaction group RG is in equilibrium by computing the fractional
        // difference of the actual and equilibrium abundances for all isotopic species
        // in the reaction group.
    
if(t < 1.04e-20) printf("\n\nt=%5.3e",t);
        for (int i = 0; i < niso; i++) {
            
            // Compute absolute value of deviation of abundances from 
            // equilibrium values for this reaction group.
            
            if(isoY[i] > 0 && isoYeq[i] > 0){
                eqcheck[i] = abs( isoY[i] - isoYeq[i] ) / max(isoYeq[i], GZ);
            } else {
                eqcheck[i] = 1e24;  // Dummy large number
            }
            eqRatio[i] = eqcheck[i]/equiTol;
           
            
            // Store some min and max values for this RG
            
            if (eqRatio[i] < minRatio) minRatio = eqRatio[i];
            if (eqRatio[i] > maxRatio) maxRatio = eqRatio[i];
            if (isoYeq[i] < Yminner) Yminner = isoYeq[i];
if(t < 1.04e-20) printf("\nRG=%d isEquil=%d i=%d %s isoY[i]=%5.3e isoYeq[i]=%5.3e eqcheck[i]=%5.3e maxRatio=%5.3e",
    RGn, isEquil, i, isolabel[i], isoY[i], isoYeq[i], eqcheck[i], maxRatio
);
            
        }
            
        // Set isEquil to false if any eqcheck[] greater than equiTol or if the 
        // time is before the time to allow equilibration equilTime, and true 
        // otherwise. useDevious and useEquilY control whether value of maxDevious or
        // value of max deviation from equil abundances (or both) determine whether
        // a RG is in equilibrium.
        
        if(useDevious && !useEquilY){
            
            if (t > equilTime && thisDevious < deviousMax) { // thisDevious condition
                addToEquilibrium();
            } else {
                isEquil = false;
            }
            
        } else if (useEquilY && !useDevious){
            
            if (t > equilTime && maxRatio < 1) {  // population condition
                addToEquilibrium();
            } else {
                isEquil = false;
            }
            
        } else if (useEquilY && useDevious){
            
            if (t > equilTime && (maxRatio < 1 || thisDevious < deviousMax)) {  // dual condition
                addToEquilibrium();
            } else {
                isEquil = false;
            }
        }
        
        // Set the activity array for each reaction in reaction group to true 
        // if not in equil and false if it is, if we are imposing equilibrium.
        
        if (doPE && t > equilTime) {
            
            for (int i = 0; i < numberMemberReactions; i++) {
                int ck = memberReactions[i];
                reacIsActive[ck] = !isEquil;
            }
        }
            
    }    // End function ReactionGroup::computeEqRatios()
    
    
    // -----------------------------------------------------------
    // Function ReactionGroup::addToEquilibrium() to add 
    // reaction group to equilibrium if it was not in equilibrium
    // in last timestep but now satisfies the equilibrium 
    // conditions for this timestep.
    // -----------------------------------------------------------
    
    void addToEquilibrium(){
        
        isEquil = true;
        totalEquilRG ++;
        
        // ******************************************************************************            
            // NOTE: Following for-loop causes Asy to produce erroneous wiggles
            // Not sure why, but comment out for now.  Maybe because the correct
            // implementation is the last loop above in computeEqRatios()?
                
                //for (int i = 0; i < numberMemberReactions; i++) {
                //    int ck = memberReactions[i];
                //    reacIsActive[ck] = false;         
                //}
        // ******************************************************************************
        
        if(showAddRemove) fprintf(pFileD,
        "\n*** ADD RG %d Steps=%d RGeq=%d t=%6.4e logt=%6.4e devious=%6.4e Rmin=%6.4f Rmax=%6.4f", RGn, totalTimeSteps, totalEquilRG, t, log10(t), thisDevious,
        minRatio, maxRatio);

    }
    
    
    // -----------------------------------------------------------
    // Function ReactionGroup::removeFromEquilibrium() to remove 
    // reaction group from equilibrium if it was in equilibrium
    // in last timestep but no longer satisfies the equilibrium 
    // conditions for this timestep.
    // -----------------------------------------------------------
    
    void removeFromEquilibrium(int where) {
        
        isEquil = false;
        //thisDevious = abs((equilRatio - kratio) / max(kratio, GZ));
        
        totalEquilRG -- ;
        
        for (int i = 0; i < numberMemberReactions; i++) {
            int ck = memberReactions[i];
            reacIsActive[ck] = true;         
        }
        
        if(showAddRemove)
        fprintf(pFileD,
        "\n*** REMOVE RG %d where=%d Steps=%d RGeq=%d t=%6.4e logt=%6.4e devious=%6.4e Rmin=%6.4f Rmax=%6.4f",
        RGn,where,totalTimeSteps,totalEquilRG,t,log10(t),thisDevious,minRatio,maxRatio);
        
    }
    
    
    // ----------------------------------------------------------------
    // Function ReactionGroup::speciesIsInRG() to determine if given 
    // species with index speciesIndex is in any of the reactions of a 
    // reaction group,where speciesIndex is the array index i for 
    // isotope quantitites like Z[i] or Y[i]. Returns true (1) if it
    // is and false (0) if not.
    // ----------------------------------------------------------------
    
    bool speciesIsInRG(int speciesIndex) { 
        
        int sindex = speciesIndex;
        
        // Loop over member reactions in the RG
        
        for(int i=0; i<numberMemberReactions; i++){

            // Loop over isotopes in reactants
            
            for (int j=0; j<numberReactants[i]; j++){ 
                if(isoindex[j] == sindex){
                    return true;
                }
            }
            
            // Loop over isotopes in products
            
            for (int j=0; j<numberProducts[i]; j++){
                if(isoindex[j+numberReactants[i]] == sindex){
                    return true;
                }
            }
        }
        
        return false;
        
    }       // End function speciesIsInRG(int)
    
    
};      // End class ReactionGroup




// ----------------------------------------------------------------
// Class Integrate with functions to integrate the reaction network.  
// Inherits from class Utilities.
// ----------------------------------------------------------------


class Integrate: public Utilities {
    
    // Make data fields private,with external access to them through public setter 
    // and getter functions.  Its static functions can be called directly from the class
    // without having to instantiate.
    
    private:
        
    public:
        
        // Static function Integrate::doIntegrationStep() that may be invokied
        // to execute a single integration step.  Assumes that all fluxes have
        // been calculated, and relevant fluxes set to zero if PE approximation
        // and the reaction group has been judged to be in equilibrium.
        // Declared static so that it can be invoked directly from the class as
        // Integrate::doIntegrationStep(), without having to instantiate an 
        // Integrate object.
        
        static void doIntegrationStep(){
            
            double sumXhalf;
            double sumXfull;
            double dE_half1;
            double dE_half2;
            double E_half1 = 0.0;
            double E_half2 = 0.0;
            double sumhalves;
            
            // Increment integration step counter for the timestep we are
            // about to execute. Time at the end of this timestep will
            // be set near the end of this function.
            
            totalTimeSteps ++; 
            
            // Choose massTol parameter
            
            double eqCut = 0.15;
            
            //if(t > 1e4){
            if(doPE && eqFrac > eqCut){
                massTol = massTol_asyPE;  // If enough equilibrated RG
            } else {
                massTol = massTol_asy;    // If too few equilibrated RG
            }

            // Store quantities from previous timestep
            
            dtLast = dt;
            sumX = Utilities::sumMassFractions();
            sumXlast = sumX;
            t_saved = t;
            dt_saved = dt;       // Will hold final dt in this timestep.
            storeCurrentY();     // For later restoration
            
            // Compute max stable forward Euler timestep
            
            if(fastestCurrentRate > 0){
                dt_FE = sf/fastestCurrentRate;
            } else {
                dt_FE = 0.00001*t;
            }
                
            
            // Compute timestep for explicit asymptotic method.  The 
            // diagnostic choice2 indicates whether dt_FE or dt_EA 
            // is chosen as the timestep.
            
            choice2 =  0;
            
            dt_EA = computeTimeStep_EA(dtLast);
            if(dt_EA < dt_FE) choice2 = 1;
            dt = min(dt_FE,dt_EA); 
            
            if(dt == 0){
                printf ("\n\n*** STOP: dt=0 in doIntegrationStep(); step=%d t=%7.5e ***\n\n",
                    totalTimeSteps,t);
                exit(1);
            }
            
            // We will estimate error by computing difference in sumX 
            // between full timestep and at end of two half timesteps. 
            
            dt_saved = dt;         // Save full dt for this step
            t = t_saved;           // Save initial time for this step
            t_end = t_saved + dt;  // This will be end time for this step
            dt_half = 0.5*dt;      // Half of full timestep
            
            // Now execute a full timestep and store sumX. Note that
            // updatePopulations(dt) updates sumX. The diagnostic 
            // flag dtMode labels whether we are computing before the
            // first integration timestep (dtMode=-1),taking the full 
            // timestep (dtMode=0),the first half timestep (dtMode = 1),
            // or the second half timestep (dtMode = 2).
            
            dtMode = 0;
            
            updatePopulations(dt);
            sumXfull = sumX;
            
            // Store the values of Y for the full step for later
            // diagnostics
            
            for(int i=0; i<ISOTOPES; i++){
                YfullStep[i] = Y[i];
            }
            
            // Now execute the first of two half-steps
            
            restoreCurrentY();
            dtMode = 1;
            updatePopulations(dt_half);
            
            // Compute energy release for first half-step
            
            dE_half1 = dE_halfstep()*dt_half;
            
            // Now update time for second half-step
            
            t = t_saved + dt_half; 

            // Recompute fluxes since Ys have changed in 1st half step
            
            computeReactionFluxes();
            dtMode = 2;
            
            // Execute second half-step.  First we
            // update Y0 array with results from first half-step
            // as new starting point for second half-step.
            
            updateY0();
            
            // Now update populations for 2nd half step. Since
            // two half steps should be more accurate than one
            // full step,we will accept this population update
            // as the final result for this integration step.
            
            updatePopulations(dt_half);
            sumXhalf = sumX; 
            
            // Compute energy release in second half step
            
            dE_half2 = dE_halfstep()*dt_half;
            
            // Set time to end of 2nd half step since we are
            // updating populations corresponding to the time
            // at the end of this integration step.
            
            t = t_end; 
            
            // Compute total energy release for both half steps based on 
            // Q-values. This energy release will be accepted as the result for this
            // integration step.
            
            sumhalves = dE_half1 + dE_half2;
            ERelease += sumhalves;
            netdERelease = (sumhalves/2.0) / dt_half;
            
            // Estimate the local error associated with the size of the
            // timestep by comparing sumX for full timestep and
            // for two successive half-steps
            
            Error_Observed = abs(sumXhalf - sumXfull);
            Error_Desired = EpsA;       // Neglect EpsR for now
            
            // Get new trial timestep for next integration step by using
            // the local error observed for this timestep compared with the
            // error desired to predict a timestep having near the local
            // error desired.

            dt = computeNextTimeStep(Error_Observed, Error_Desired, dt_saved);
            dtMode = -1;    // Indicates not in the integration step
            
        }    // End of doIntegrationStep
        
        
        
        // Function Integrate::storeCurrentY() to store current Ys. Restore 
        // using restoreCurrentY().
        
        static void storeCurrentY(){
            for(int i=0; i<ISOTOPES; i++){
                Ystore[i] = Y[i];
            }
        }
        
        // Function Integrate::restoreCurrentY()to restore Ys that were stored 
        // using storeCurrentY()
        
        static void restoreCurrentY(){
            for(int i=0; i<ISOTOPES; i++){
                Y[i] = Ystore[i];
            }
        }
        
        // Function Integrate::computeTimeStep_EA (double,double)
        // to update explicit asymptotic timestep

        static double computeTimeStep_EA(double dt0){
    
            double dtt = dt0;
            double dtt_0;

            updatePopulations(dtt);
            
            // Iterate timestep downward if necessary to satisfy the
            // particle number conservation condition
            
            iterations = 0;
            
            diffXzero = diffX;
            
            while( diffX > massTol && iterations <= maxit){
                
                dtt_0 = dtt;
                dtt *= downbumper;
                dt = dtt;
                updatePopulations(dtt);
                diffX = abs(sumX-sumXlast);
                iterations ++;
                totalIterations ++;
                
                if(iterations == maxit){
                    printf("\n\n*** HALT: t=%6.3es dt_iterate=%6.3es; iterations=maxit=%d", 
                        t, dt, iterations);
                    if(fixCNO) printf("\n          fixingCNO_now=%d startX_fixCNO_time=%6.3es", 
                        fixingCNO_now, startX_fixCNO_time);
                    printf("\n\n");
                    exit(1);
                }
            }
            
            if (iterations > mostIterationsPerStep){
                mostIterationsPerStep = iterations;
                maxIterationStep = totalTimeSteps;
                maxIterationTime = t;
            }
            
            diffXfinal = diffX;
            
            // Ensure timestep not larger than fraction dt_ceiling of time,
            // for accuracy. 
            
            dtt = min(dtt, dt_ceiling*t);
            
            // If timestep would cause t+dt to be larger than the next
            // plot output step, reduce trial dt to be equal to upfac times
            // the next plot output step.
            
            
            
            dt_desired = dtt;
            double upfac = 1.5;
            double gap = upfac*nextPlotTime - t_saved;
            
            if(dtt > gap && gap > 0){
                dtt = gap;
            }
            
//             dt_desired = dtt;
//             double upfac = 0.0001;
//             double gap = nextPlotTime - t_saved;
//             
//             if(dtt > gap && gap > 0){
//                 dtt = (1+upfac)*gap;
//             } 
            
//             dt_desired = dtt;
//             double upfac = 1.5;
//             double gap = nextPlotTime - t_saved;
// 
//             if(dtt > gap && gap > 0){
//                 dtt = upfac*gap;
//             } 
            
            // dtt should be a positive number.  Exit if it isn't.
            
            if(dtt <= 0){
                printf("\n\n*** STOP: dt=0 in computeTimeStep_EA(); step=%d,t=%7.5e\n\n",
                    totalTimeSteps,t);
                exit(1);
            }
            
            return dtt;
            
        }  // End of computeTimeStep_EA
        
        
        // Function Integrate::computeNextTimeStep(double,double,double)
        // to compute the next timestep based on the error observed in the
        // previous timestep.
        
        static double computeNextTimeStep(double error, double error_D, double dt_old){
            
            E_R = error/error_D; 
            
            if(E_R > 0.1){
                dt_new =  0.9 * dt_old / E_R;
                choice1 = 0;
            } else if(E_R < 0.5) {
                dt_new =  0.5 * dt_old / sqrt(E_R);   // Note: Divides by 0 if E_R = 0
                choice1 = 1;
            } else {
                dt_new = dt_old;
                choice1 = 2;
            }
            
            // Restrict dt_new to not be larger than twice old timestep
            // and not smaller than half old timestep
            
            double maxupdt = 2.0;
            double maxdowndt = 0.5;
            double dt_MAXFAC = 0.1;
            
            dtmin = min(dt_new, maxupdt*dt_old);
            dt_new = max( dtmin, maxdowndt*dt_old );
            
            // Don't let dt exceed dt_MAXFAC * t for accuracy reasons
            
            dt_new = min(dt_new, dt_MAXFAC*t);
            
            // Return the new trial timestep that will be the starting point
            // for the next integration step.
            
            return dt_new;
        }
        
        
        // Function Integrate::updatePopulations(double) to do population update 
        // with current dt and current fluxes.
        
        static void updatePopulations(double dtt){

            // If using the QSS approximation, apply the QSS approximation 
            // to all isotopes
            
            if(doQSS){
                QSSupdate();
            }
            
            // If using asymptotic approximation,determine which species satisfy the
            // asymptotic condition. 
            
            if(doASY){
                
                // Determine which isotopes satisfy asymptotic condition
                
                int numAsy = 0;
                
                for(int i=0; i<ISOTOPES; i++){
                    isAsy[i] = checkAsy(keff[i]);
                }
                
                dt = dtt;
                
                // If Asy or Asy+PE,compute with Asy+FE algorithm (with fluxes removed
                // by PE approximation (individual fluxes in RG that are in equilibrium)
                // if Asy+PE.
                
                updateAsyEuler();
            }
            
//             if(fixCNO && X[index1H] < startX_fixCNO){
//                 correctCNOCycle();
//             }
            
        }  // End of updatePopulationsa
        

        // The function Integrate::updateAsyEuler() uses the fluxes
        // to update the populations for this timestep. We determine whether each isotope 
        // satisfies the asymptotic condition. If it does we update with the asymptotic formula. 
        // If not, we update numerically using the forward Euler formula. 
        
        static void updateAsyEuler(){
            
            sumX = 0.0;
            
            for(int i=0; i<numberSpecies; i++){	
                
                if(isAsy[i]){
                    Y[i] = asymptoticUpdate(i,FplusSum[i],FminusSum[i],Y0[i],dt);
                } else {
                    Y[i] = eulerUpdate(i,FplusSum[i],FminusSum[i],Y0[i],dt);
                }
                X[i] = Y[i] * (double) AA[i];
                sumX += X[i];
                
            }
            
            diffX = abs(sumX-sumXlast);
            
        }    // End function updateAsyEuler()
    
     
     
     // Function Integrate::eulerUpdate(int, double, double, double, double) to update 
     // by the forward Euler method. Returns the updated value of Y.
        
    static double eulerUpdate(int i, double fplus, double fminus, double y0, double dtt){
        
        double newY;
        double dY;
        
        dF[i] = fplus - fminus;
        dY = dF[i]*dtt;
        newY = y0 + dY;
        return newY;            // New Y for forward Euler method
        
    }
    
    // Function Integrate::asymptoticUpdate(int double, double, double, double) to update
    // by the asymptotic method using Sophia He formula. Returns the updated value of Y.
    
    static double asymptoticUpdate(int i, double fplus, double fminus, double y, double dtt){
        
        double newY = (y + fplus*dtt)/(1.0 + fminus*dtt/y);  
        return newY;
        
    }
    
    // Function Integrate::checkAsy() to determine whether an isotope satisfies the
    // asymptotic condition. Returns true (1) if it does and false (0) if not.
    // Two overloaded forms taking either Fminus and Y0,or keff as arguments.
    
    static bool checkAsy(double Fminus,double YY){
        
        asycheck = Fminus*dt/YY;
        if(YY > 0.0 && asycheck > 1.0){
            return true;
        } else {
            return false;
        }
        
    }
    
    static bool checkAsy(double k){
        
        asycheck = k*dt;
        if(asycheck > 1.0){
            return true;
        } else {
            return false;
        }
        
    }
    
    
    // Function to update by the Quasi-Steady-State (QSS) approximation. This function 
    // replicates method steadyState(dt) in Java code. nitQSS is the number of
    // predictor-corrector iterations.  Found with Java code that increasing beyond 1
    // did not help much,so generally use nitQSS = 1,but code allows nitQSS > 1.
    
    static void QSSupdate(){
        
        // Iteration loop (nitQSS is number of predictor-corrector iterations; 
        // normally nitQSS=1)
        
        int nitQSS = 1;
        
        for (int i = 0; i < nitQSS; i++) {
            
            ssPredictor();
            ssCorrector();
            
            // Recompute fluxes if another predictor-corrector iteration follows
            
            if (nitQSS > 1) {
                computeReactionFluxes();
            }
        }
    }
    
    
    // ----------------------------------------------------------------------
    // Integrate::ssPredictor() to implement steady-state predictor step
    // ----------------------------------------------------------------------
    
    static void ssPredictor() {
        
        // Save current values of F+,F- and keff for later use. Y0[]
        // already contains the saved populations before update 
        // (i.e.,from last timestep)
        
        for (int i = 0; i < ISOTOPES; i++) {

            FplusZero[i] = FplusSum[i];
            keffZero[i] = keff[i];
            
        }
        
        // Loop over all active isotopes and calculate the predictor
        // populations. Unlike for asymptotic method where we update with
        // the asymptotic approximation if kdt >= 1 and with forward euler if
        // if kdt < 1,we will update isotope abundances with the same 
        // QSS predictor-corrector,irrespective of value of keff*dt for 
        // an isotope.
        
        for (int i = 0; i < ISOTOPES; i++) {
            
                double kdt = keff[i] * dt;
                double deno = 1.0 + kdt*alphaValue(kdt);
                
                Y[i] = Y0[i] + (FplusSum[i] - FminusSum[i])*dt / deno;
                X[i] = Y[i] * (double)AA[i];

        }
        
        // Update all fluxes for the corrector step
        
        computeReactionFluxes();
        
    }
    
    
    // ------------------------------------------------------------------
    // Integrate::ssCorrector() to implement steady-state corrector step
    // based on earlier predictor step
    // ------------------------------------------------------------------
    
    static void ssCorrector() {
        
        double kBar;
        double kdt;
        double alphaBar;
        double FplusTilde;
        sumX = 0.0;
        
        // Loop over all isotopes and update the result of the predictor step
        
        for (int i = 0; i< ISOTOPES; i++) {

            kBar = 0.5 * (keffZero[i] + keff[i]);
            kdt = kBar * dt;
            alphaBar = alphaValue(kdt);
            FplusTilde = alphaBar * FplusSum[i] + (1.0 - alphaBar) * FplusZero[i];
            Y[i] = Y0[i] + ((FplusTilde - kBar * Y0[i]) * dt) / (1 + alphaBar * kdt);
            X[i] = Y[i] * (double)AA[i];
            
            // For reference,keep track of the isotopes that would satisfy the
            // asymptotic condition if we were using the asymptotic approximation
            // instead of QSS
            
            if (kdt >= 1.0) {
                
                isAsy[i] = true;
                totalAsy ++;
                
            } else {
                
                isAsy[i] = false;
                
            }
            
            // Update sum of mass fractions
            
            sumX += X[i];
                
        }
    }
 
    
    // ----------------------------------------------------------------
    // Function Integrate::alphaValue(double) to calculate alpha(kdt) 
    // for steady state approximation.
    // ----------------------------------------------------------------
    
    static double alphaValue(double a) {
        
        double aa = a;
        
        // Following necessary to start integration correctly.
        
        if (aa < 1.e-20) aa = 1e-20; 
        
        double ainv = 1.0/aa;
        double a2 = ainv * ainv;
        double a3 = a2 * ainv;
        
        return (180.0 * a3 + 60.0 * a2 + 11.0 * ainv + 1.0)
            / (360.0 * a3 + 60.0 * a2 + 12.0 * ainv + 1.0);
        
	}
    
    
};    // End class Integrate


//----------------END CLASS DEFINITIONS ----------------




// Declare pointer used to access the fields and functions of class Species

Species *SpeciesPtr;

// Create an array of Species objects isotope[] to hold information and functions
// for the isotopic species in the network. Each element of the array will be
// a Species object corresponding to a different isotope of the network.

Species isotope[ISOTOPES];

// Declare pointer used to access the fields and functions of class Reaction

Reaction *ReactionPtr;

// Create an array of Reaction objects reaction[] to hold information and functions
// for the reactions in the network. Each element of the array will be
// a Reaction object corresponding to a different reaction of the network.

Reaction reaction [SIZE];

// Create pointer to an array of ReactionGroup objects RG[] to hold information and 
// functions for reactions groups in the network. Each element of the array will be
// a ReactionGroup object corresponding to a different reaction group of the network.
// Memory for array will be allocated dynamically below to the size
// given by numberRG (which is computed in the class ReactionVector)

ReactionGroup *RG;   // Pointer to 1D array for reaction groups



// ---------------------------------
// ------- Main CPU routine --------
// ---------------------------------


int main() { 

    // Open file to output network info
    
    pfnet = fopen("gnu_out/network.data","w");
    
    // Open a file for diagnostics output
    
    pFileD = fopen("gnu_out/diagnostics.data","w");
    
    // Write the time
    
    fprintf(pFileD, Utilities::showTime());
    printf("\n%s", Utilities::showTime());
    
    // Initialize the array of current partition functions to unity.
    
    fill_n (currentPF, ISOTOPES, 1.0);
    
    // Set labels and check consistency of choice for explicit algebraic methods set.
    // Generally we use either asymptotic (Asy) or quasi-steady-state (QSS) algorithms.
    // In either case we may choose to add the partial equilibrium (PE) algorithm. So
    // valid options are Asy, QSS, Asy+PE, and QSS+PE.
    
    if(doASY && !doPE){
        doQSS = false;
        intMethod = "ASY method";
        fprintf(pFileD,"Using ASY method\n");
    } else if (doQSS && !doPE) {
        intMethod = "QSS method";
        doASY = false;
        fprintf(pFileD,"Using QSS method\n");
    } else if (doASY && doPE){
        intMethod = "ASY+PE method";
        fprintf(pFileD,"Using ASY+PE method\n");
    } else if (doQSS && doPE){
        intMethod = "QSS+PE method";
        fprintf(pFileD,"Using QSS+PE method\n");
    }

    // Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
    // realistic calculation the temperature and density will be passed from the hydro 
    // code in an operator-split coupling of this network to hydro. Here we hardwire
    // them for testing purposes.  These will be used to calculate the reaction
    // rates in the network. If we assume operator splitting, the temperature 
    // and density are assumed constant for each network integration. We also allow
    // the possibility below to interpolate the temperature and density from a
    // hydrodynamical profile as a function of time.
    
    T9 = T9_start;
    rho = rho_start;
    
    // Display basic parameters
    
    showParameters();
    
    // Initialize reacIsActive[] array to true and isEquil field of reaction[]
    // objects to false.
    
    for (int i=0; i<SIZE; i++){ 
        reacIsActive[i] = true;
        reaction[i].setisEquil(false);
        //reaction[i].setreacGroupSymbol(Utilities::stringToChar(RGstring[i]));
    }
    
    // Determine whether the network contains Be-8, which is handled as
    // two alpha particles because of the rapid decay compared to network 
    // timescales.
    
    for(int i=0; i<ISOTOPES; i++){
        
        hasBe8 = false;
        hasAlpha = false;
        
        if(isotope[i].getZ()==4 && isotope[i].getN()==4){
            hasBe8 == true;
            indexBe8 == i;
        }
        
        if(isotope[i].getZ()==2 && isotope[i].getN()==2){
            hasAlpha == true;
            indexAlpha == i;
        }
    }
        
    
    // Read in network file and associated partition functions.  This is required only
    // once at the beginning of the entire calculation.  
    
    char *networkFilePtr = networkFile;
    readNetwork(networkFilePtr);
    writeNetwork();
    
    // Read in rate library data from a file. This is required only once, at the
    // beginning of the entire calculation.
    
    char *rateLibraryFilePtr = rateLibraryFile;
    readLibraryParams(rateLibraryFilePtr);
    
    // If using a hydrodynamical profile,read in the file containing
    // the hydro profile and store variables.
    
    if(hydroProfile){
        char *hydroFilePtr = hydroFile;
        readhydroProfile(hydroFilePtr);
    }
    
    // Print out some quantitites from the Reaction object reaction[].  
    
    for(int i=0; i<SIZE; i++){
        
        fprintf(pfnet,
            "\n%d %s reacClass=%d #react=%d #prod=%d isEC=%d isReverse=%d Q=%5.4f prefac=%5.4f RGsymb:%s",
            reaction[i].getreacIndex(),
            reacLabel[i], 
            reaction[i].getreacClass(),
            reaction[i].getnumberReactants(),
            reaction[i].getnumberProducts(),
            reaction[i].getisEC(),
            reaction[i].getisReverse(),
            reaction[i].getQ(),
            reaction[i].getprefac(),
            Utilities::stringToChar(RGstring[i])
            //reaction[i].getreacGroupSymbol()   // Always returns a<->b?
        );
    }

    // Find the time intervals for plot output during the integration. After this
    // function is executed the plotSteps target time intervals for output will
    // be in the array plotTimeTargets[]. In the integration the ith output step will 
    // be triggered as soon as the time t is >= plottimeTargets[i].  The actual time of the output
    // (which will usually be slightly larger than plottimeTargets[i]) will be stored in tplot[i].
    // The variables start_time and stop_time define the range of integration.  The variable
    // startplot_time allows the plotting interval output in to be a subset of
    // the full integration interval.
    
    Utilities::log10Spacing(max(start_time,startplot_time),stop_time,plotSteps,plotTimeTargets);
    
    // Find for each isotope all reactions that change its population.  This analysis of
    // the network is required only once at the very beginning of the calculation (provided
    // that the network species and reactions remain the same for the entire calculation).
    // The work is done by the function ReactionVector::parseF().  First allocate array
    // memory.
    
    // Number of F+ and F- components for each isotope
    
    numFluxPlus = (int*) malloc(sizeof(int) * numberSpecies);
    numFluxMinus = (int*) malloc(sizeof(int) * numberSpecies);
    
    // Arrays for temporary storage
    
    tempInt1 = (int*) malloc(sizeof(int) * numberSpecies * numberReactions/2);
    tempInt2 = (int*) malloc(sizeof(int) * numberSpecies * numberReactions/2);
    
    /*
     * Use ReactionVector::parseF() to find the contributions to F+ and F- of each reaction for each 
     * isotope. This is executed only once at the beginning of the entire calculation to determine 
     * the structure of the network.  The function is static so it can be called directly from the
     * class without instantiating.
     */
    
    ReactionVector::parseF();
    
    // Print out the network species vector
    
    if(showDetails) fprintf(pFileD,
    "\n\nNETWORK SPECIES VECTOR (%d components):\n\nIndex  Species    Z     N", numberSpecies);
    
    for(int i=0; i<numberSpecies; i++){
        if(showDetails) fprintf(pFileD,
            "\n%5d    %5s  %3d  %4d", i, isoLabel[i], Z[i], N[i]);
    }
    
    if(showDetails) fprintf(pFileD,"\n");
    
    // Use the information gleaned from ReactionVector::parseF() to define 
    // the reaction vectors for the network using the static makeReactionVectors 
    // function of the class ReactionVector.
    
    printf("\n\nMaking Reaction Vectors ...");
    cout.flush();
    
    ReactionVector::makeReactionVectors();
    
    // Use static function ReactionVector::sortReactionGroups() to sort 
    // reactions into partial equilibrium reaction groups by comparing 
    // reaction vectors.
    
    printf("\nSorting Reaction Vectors ...");
    cout.flush();
    
    ReactionVector::sortReactionGroups();
    
    // Allocate dynamically memory for an array of ReactionGroup objects 
    // of dimension numberRG, where numberRG was determined by 
    // ReactionVector::sortReactionGroups() above.
    
    printf("\nAllocating Reaction Group Objects ...");
    cout.flush();
    
    RG = (ReactionGroup*) malloc(sizeof(ReactionGroup) * numberRG);
    
    // Create ReactionGroup objects RG[] and assign values for various fields
    // using the function assignRG().
    
    assignRG();
    
    // Populate boolean array RGisoMembers[RG][isotopes] giving the isotopes 
    // appearing in each RG. Entry in array is true (1) if the isotope appears
    // in the corresponding reaction group and false (0) otherwise.
    
    for(int i=0; i<numberRG; i++){
        
        for(int j=0; j<ISOTOPES; j++){
            if( RG[i].speciesIsInRG(j) ){
                RGisoMembers[i][j] = true;
            } else {
                RGisoMembers[i][j] = false;
            }
        }
        
    }
    
    // Allocate memory for 1D arrays to hold non-zero F+ and F- for all reactions 
    // for all isotopes, the arrays holding the species factors FplusFac and 
    // FminusFac,and also arrays to hold their sums for each isotope. 
    // ReactionVector::parseF() must be run first because it determines 
    // totalFplus and totalFminus.
    
    Fplus = (double*) malloc(sizeof(double) * totalFplus);
    Fminus = (double*) malloc(sizeof(double) * totalFminus);
    FplusFac = (double*) malloc(sizeof(double) *totalFplus);
    FminusFac = (double*) malloc(sizeof(double) * totalFminus);
    FplusSum = (double*) malloc(sizeof(double) * numberSpecies);
    FminusSum = (double*) malloc(sizeof(double) * numberSpecies);
    dF =  (double*) malloc(sizeof(double) * numberSpecies);
    
    // Allocate memory for arrays that hold the index of the boundary between 
    // different isotopes in the Fplus and Fminus 1D arrays. 
    
    FplusMax = (int*) malloc(sizeof(int) * numberSpecies);
    FplusMin = (int*) malloc(sizeof(int) * numberSpecies);
    FminusMax = (int*) malloc(sizeof(int) * numberSpecies);
    FminusMin = (int*) malloc(sizeof(int) * numberSpecies);

    // Allocate memory for 1D arrays that will be used to map finite F+ and F- to 
    // the Flux array.
    
    FplusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
    FminusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
    FplusIsotopeIndex = (int*) malloc(sizeof(int) * totalFplus);
    FminusIsotopeIndex = (int*) malloc(sizeof(int) * totalFminus);
    
    // Allocate memory for 1D arrays that will hold the index of the isotope for 
    // the F+ or F- term
    
    MapFplus = (int*) malloc(sizeof(int) * totalFplus);
    MapFminus = (int*) malloc(sizeof(int) * totalFminus);
    
    // Call function Reaction::setupFplusFminus() to set up F+ and F- index for each
    // isotope and to find non-vanishing F+ and F- source terms in network.
    // Function declared static so it can be called as Reaction::setupFplusFminus() 
    // without having to instantiate.
    
    Reaction::setupFplusFminus();
    
    
    
    // -----------------------------------------------------
    // *** Set up for main time-integration while-loop ***
    // -----------------------------------------------------
    
    // Check if IS0TOPES and numberSpecies not equal. Normally they should be
    // as present code is set up.
    
    if (ISOTOPES != numberSpecies){
        printf("\n\n***** CAUTION:  ISOTOPES=%d numberSpecies=%d *****\n\n",
            ISOTOPES,numberSpecies);
    }
    
    // Treat CNO closed cycle
    
    if(CNOinNetwork = checkForCNOCycle()){
        
        //printf("\n\n+++++We have CNO CNOinNetwork=%d", CNOinNetwork);
        
        indexCNOCycle();
        
        int isosymb = Utilities::returnNetIndexSymbol(isoLabel[1]);
        //printf("\nreturnNetIndexSymbol: index=%d %s", isosymb, isoLabel[1]);
        
        int reacy = Utilities::returnNetIndexSymbol(*reacLabel+1);
        
        string sreacy = Utilities::charArrayToString(reacLabel[1], LABELSIZE);
        
        //printf("\n%s %s", Utilities::stringToChar(sreacy), reacLabel[1]);
        
        Utilities::compareTwoCharArrays(Utilities::stringToChar("p+n15-->he4+c12"), reacLabel[21]);
        
        
    }
    
    
    
    printf("\n\nBEGIN TIME INTEGRATION:\n");
    fprintf(pFileD,"\n\n\n\n                 --- BEGIN TIME INTEGRATION ---\n");
    
    dt = dt_start;              // Integration start time
    t = start_time;             // Current integration time
    totalTimeSteps = 1;         // Integration step counter
    totalEquilRG = 0;           // Number quilibrated reaction groups
    totalEquilReactions = 0;    // Number equilibrated reactions
    totalAsy = 0;               // Number asymptotic species
    ERelease = 0.0;             // Total E released from Q values
    plotCounter = 1;            // Plot output counter
    fastestOverallRate = 0.0;   // Initialize fastest overall rate
    timeMaxRate = 0.0;          // Initialize slowest overall rate
    
    lastNetworkMass = returnNetworkMass();  // Initialize total mass of network
    
    // Compute initial rates if hydroProfile = false. Rates won't
    // change in the integration and don't need to be computed again.  If either
    // T9 or rho change, the rates will be recomputed at each integration step.
    // Use functions of Reaction class to compute reaction rates. We have instantiated
    // a set of Reaction objects in the array reaction[i], one entry for each
    // reaction in the network. Loop over this array and call the computeRate()
    // function of Reaction on each object. 
    
    if(!hydroProfile){
        
        updatePF();   // Update partition functions for initial T9 if constant T
        
        for(int i=0; i<SIZE; i++){
            reaction[i].computeConstantFacs(T9,rho);
            reaction[i].computeRate(T9,rho);
        }
    }
    
    // Summarize computed rates
    
    if(!hydroProfile){
        
        if(showDetails) fprintf(pFileD,"\nINITIAL COMPUTED RATES\n");
        for(int i=0; i<SIZE; i++){
            reaction[i].showRates();
        }
    }
    
    if(!hydroProfile){
        if(showDetails) fprintf(pFileD,
            "\n\n**** Rates won't be computed again since T and rho won't ");
        if(showDetails) fprintf(pFileD,
            "change in integration ****\n");
    } else {
        if(showDetails) fprintf(pFileD,
            "\n\n**** Rates will be recomputed at each timestep since T and ");
        if(showDetails) fprintf(pFileD,"rho may change ****\n");
    }
    
    // Instantiate hydro temperature interpolator object
    
    SplineInterpolator interpolateT = SplineInterpolator (hydroLines, hydroTime, hydroTemp);
    SplineInterpolator interpolateRho = SplineInterpolator (hydroLines, hydroTime, hydroRho);
    
    // Initialize interpolator objects if using a hydro profile
    
    if(hydroProfile){
        
        interpolateT.spline(hydroTime, hydroTemp, hydroLines, hydroLines);
        interpolateRho.spline(hydroTime, hydroRho, hydroLines, hydroLines);
    }
    
    // Open files for plot output. Assumes that the subdirectory
    // gnu_out already exists. If it doesn't, will compile but
    // may crash when executed.
    
    plotfile1 = fopen("gnu_out/plot1.data", "w");   // t,dt,E,dE,X
    plotfile2 = fopen("gnu_out/plot2.data", "w");   // dt,Rmax
    plotfile3 = fopen("gnu_out/plot3.data", "w");   // fluxes
    plotfile4 = fopen("gnu_out/plot4.data", "w");   // hydro profile
    plotfile5 = fopen("gnu_out/plot5.data", "w");   // time & energy values, not log

    // Setup files for plot output during the integration by writing headers

    plotFileSetup();
    
    
    // ------------------------------------ //
    // *** Begin main integration loop ***  //
    // ------------------------------------ //
    
    
    totalIterations = 0;
    XcorrFac = 1.0;
    reacLibWarn = false;
    
    Utilities::startTimer();    // Start a timer for integration
    
    while(t < stop_time){ 
        
        if(plotCounter > plotSteps) {
            printf("\n\nStopping Integration: plotCounter=%d t=%7.4e dt=%7.4e",
                plotCounter, t, dt);
            break;
        }
        
        // Next target plot output time.  Use to keep chosen dt from being
        // much larger than the time to next plot output, which can occur
        // if the number of plot steps plotSteps is large (especially at
        // early times in the integration).
        
        nextPlotTime = plotTimeTargets[plotCounter-1];
        
        // Initialize fastest and slowest rates for this timestep
        
        fastestCurrentRate = 0.0;
        slowestCurrentRate = 1e30; 
        
        // Update Y0[] array with current Y[] values
        
        updateY0();
        
        // Specify temperature T9 and density rho. If hydroProfile = false, a constant
        // temperature is assumed for the entire network calculation, set by T9_start
        // above.  Otherwise (hydroProfile = true) we here interpolate the temperature
        // from a hydrodynamical profile for each timestep.  Likewise for the density.
        // Since the arrays holding the hydro profile have entries in terms of log10,
        // interpolate in log10(t) if hydroProfile = true.
        
        if(hydroProfile){
            
            logTnow = interpolateT.splint(log10(t));
            logRhoNow = interpolateRho.splint(log10(t));
            
        } else {
            
            logTnow = log10(T9_start*1e9);
            logRhoNow = log10(rho_start);
        }
        
        // But then convert also to temperature in units of 10^9 K and density 
        // in units of g/cm^3 for the calculation
        
        T9 = pow(10,logTnow)/1e9;
        rho = pow(10,logRhoNow);
        
        updatePF();    // Update partition functions for current T9
    
        // Use functions of Reaction class to compute reaction rates. We have instantiated
        // a set of Reaction objects in the array reaction[i], one entry for each
        // reaction in the network. Loop over this array and call the computeRate()
        // function of Reaction on each object. If hydroProfile is false,
        // the rates only need be computed once as they won't change over this
        // integration. Otherwise they are recomputed for each integration step in
        // the following conditional.
        
        if(hydroProfile){
            
            for(int i=0; i<SIZE; i++){
                reaction[i].computeConstantFacs(T9, rho);
                reaction[i].computeRate(T9, rho);
            }
            
        }
        
        // Compute reaction fluxes corresponding to current rates and populations.
        
        computeReactionFluxes();
        
        // Find max dY/dt and corresponding isotope for this timestep based on the initial
        // value of data for the timestep. Presently for information only.
        
        getmaxdYdt();
        
        // Perform an integration step using the static method doIntegrationStep() of
        // the class Integrate.
        
        Integrate::doIntegrationStep();
        
        // Variable t now holds the time at the end of the timestep just executed.
        
        if(fixCNO && X[index1H] < startX_fixCNO){
            correctCNOCycle();
        }
        
//         // Try removing stiffness associated with beta decay in main CNO cycle by 
//         // computing the 15N and 13C populations at equilibrium from the 14N population,
//         // by requiring that in the closed cycle the rates must be equal around the 
//         // cycle at equilibrium. The boolean fixCNO controls whether this fix
//         // is applied if X[H] < startX_fixCNO.
//         
//         // 15N population
//         
//         //double fluxCycle6 = (Rate[17] + Rate[18])/(Rate[10] + Rate[11]);    // CNO
//         double fluxCycle6 = (Rate[25] + Rate[26])/(Rate[14] + Rate[15]);  // extended CNO
//         double Ycycle6 = fluxCycle6 * Y[5];
//         double Yratio6 = Ycycle6/Y[6];
//         
//         // 13C population
//         
//         //double fluxCycle3 = (Rate[17] + Rate[18])/(Rate[12] + Rate[13]);  // CNO
//         double fluxCycle3 = (Rate[25] + Rate[26])/(Rate[18] + Rate[19]);  // extended CNO
//         double Ycycle3 = fluxCycle3 * Y[5];
//         double Yratio3 = Ycycle3/Y[3];
//         
//         if(fixCNO && X[0] < startX_fixCNO){
//             
//             if(!fixingCNO_now){
//                 startX_fixCNO_time = t;
//                 fixingCNO_now = true;
//             }
//             
//             // 15N
//             
//             Y[6] = Ycycle6;
//             X[6] = Y[6] * 15;
//             
//             // 13C
//             
//             Y[3] = Ycycle3;
//             X[3] = Y[3] * 15;
//         }
        
        // Store true sumX before any renormalization.
        
        sumXtrue = Utilities::sumMassFractions();
        
        // Convert all Be-8 to alpha particles since lifetime of Be-8 to decay
		// to two alpha particles is short compared with typical integration steps.
        
        if(hasBe8 && hasAlpha){
            restoreBe8();
        }
        
        // Compute equilibrium conditions for the state at the end of this timestep 
        // (starting time for next timestep) if partial equilibrium is being 
        // implemented (doPE = true).
        
        if( (doPE && t > equilTime)){
            
            for(int i = 0; i < numberRG; i++) {
                RG[i].computeEquilibrium();
            }
            
            if(totalEquilRG > 0){
                
                // Restore species in equilibrium to their unperturbed equilibrium values
                
                restoreEquilibriumProg();
                
            }
        }
        
        // If showPE == true (so doPE == false), count RG that would be in equilibrium 
        // if doPE were true in current Asy or QSS calculations w/o PE.
        
        if(showPE){
            int totalPseudoEquilRG = 0;
            for(int i = 0; i < numberRG; i++) {
                RG[i].computeEquilibrium();
                if(RG[i].getisEquil()) totalPseudoEquilRG ++;
            }
        }
        
        // Count total asymptotic species
        
        totalAsy = 0;
        for(int i=0; i<ISOTOPES; i++){
            if (isAsy[i]){totalAsy ++;}
        }
        
        // Compute fraction of isotopes satisfying the asymptotic condition,
        // and fraction of reaction groups satisfying equilibrium condition.
        
        asyFrac = (double)totalAsy/(double)ISOTOPES;
        eqFrac = (double)totalEquilRG/(double)numberRG;
        
        // Compute energy release from mass differences
        
        networkMassDifference();
        
        
        // ---------------------------------------------------------------------------------
        // Display and output to files updated quantities at plotSteps times corresponding
        // to (approximately) equally-spaced intervals in log_10(time). The target output 
        // times are generated by Utilities::log10Spacing() and stored in the array
        // plotTimeTargets[plotSteps]. A screen and plot file output is triggered if
        // t >= plotTimeTargets[plotCounter-1],so the actual output times may be slightly
        // larger than the target output times. Plots are made with respect to actual,not
        // target, output times.
        // ---------------------------------------------------------------------------------
        
        if(t >= plotTimeTargets[plotCounter-1]){
            
            // Record int step when plotting starts
            
            if(plotCounter == 1) {totalTimeStepsZero =  totalTimeSteps;}
            
            // Output plot data to files for this timestep
            
            toPlotNow();
            
            // Output to screen for this plot step
            
            ts = "\n%d it=%d t=%6.2e dt=%6.2e dt'=%6.2e int=%d asy=%4.2f ";
            ts += "eq=%4.2f sX=%-4.3f Xfac=%-4.3f dEA*dt=%6.2e dEA=%6.2e E=%6.2e EA=%6.2e E_R=%6.2e c1=%d";
            //ts += "eq=%4.2f sX=%-4.3f Xfac=%-4.3f dE=%6.2e dEA=%6.2e E=%6.2e EA=%6.2e E_R=%6.2e c1=%d";
            ts += " c2=%d fast=%d Q=%4.2f";
            ts += " dev=%4.2e lT=%4.3f lrho=%4.2f";
            
            printf(Utilities::stringToChar(ts),
                   plotCounter, iterations, t, dt, dt_desired, totalTimeSteps,
                   asyFrac, eqFrac, sumX, XcorrFac, ECON*dEReleaseA*dt, ECON*dEReleaseA,
                   //asyFrac, eqFrac, sumX, XcorrFac, ECON*netdERelease, ECON*dEReleaseA,
                   ECON*ERelease, ECON*EReleaseA, E_R, choice1, choice2,
                   fastestCurrentRateIndex, reaction[fastestCurrentRateIndex].getQ(),
                   mostDevious, logTnow, logRhoNow
            );
            
            // Above printf writes to a buffer and the buffer is written to the screen only
            // after the buffer fills. Following command flushes the print buffer
            // after each timestep so the output to screen is continuous as the integration
            // proceeds instead of in large chunks each time the buffer fills.  This is
            // useful for diagnostics and development, but could impair efficiency of
            // code somewhat relative to letting the code decide when to empty buffer? Of
            // course in production code we would not be writing significant output to the 
            // screen anyway.
            
            cout.flush();   // Force buffer dump 
  
            // Increment the plot output counter for next graphics output
            
            plotCounter ++;
            
        }
    
    }   // End time integration while-loop
    
    displayRGstatus();   // Send status of RG at end to network.data
    
    // Close plot output files
    
    fclose(plotfile1);
    fclose(plotfile2);
    fclose(plotfile3);
    fclose(plotfile4);
    fclose(plotfile5);
    
    
    // Write parameters at end of integration. Also stops timer started
    // at beginning of the integration while-loop.
    
    showParameters();
   
    // Display final abundances and mass fractions

    printf("\nFINAL ABUNDANCES Y AND MASS FRACTIONS X\n");
    fprintf(pFileD,"\n\nFINAL ABUNDANCES Y AND MASS FRACTIONS X\n");

    for(int i=0; i<ISOTOPES; i++){
        string asylabel = "(Asy)";
        if(!isAsy[i]) asylabel = " ";
        printf("\n%d %s Y=%7.4e X=%7.4e %s",
            i, isotope[i].getLabel(), Y[i], X[i], Utilities::stringToChar(asylabel)
        );
        fprintf(pFileD,"\n%d %s Y=%7.3e X=%7.3e %s",
            i, isotope[i].getLabel(), Y[i], X[i], Utilities::stringToChar(asylabel)
        );
    }

    printf("\n\n");
    
    // Output the hydro profile read in to gnu_out/hydroProfileInput.data. 
    // (The interpolated hydro profile is output to gnu_out/plot4.data.)
    
    if(hydroProfile && plotHydroProfile) Utilities::outputHydroProfile();
    
    // Flush data buffers to force immediate output
    
    cout.flush();
    
    
    
    // **************************************************
    // To test various function calls, insert code from 
    // functionTests.cpp here
    // **************************************************

    // Close output diagnostic files
    
    fclose (pFileD);
    fclose (pfnet);
   
    // Free allocated memory

    free(Fplus);
    free(Fminus);
    free(FplusFac);
    free(FminusFac);
    free(FplusSum);
    free(FminusSum);
    free(dF);
    free(FplusMin);
    free(FplusMax);
    free(FminusMin);
    free(FminusMax);
    free(FplusIsotopeCut);
    free(FminusIsotopeCut);
    free(MapFplus);
    free(MapFminus);
    free(numFluxPlus);
    free(numFluxMinus);
    free(tempInt1);
    free(tempInt2);
    free(FplusIsotopeIndex);
    free(FminusIsotopeIndex);
    free(RG);
    
    gsl_vector_free(abundances);
    gsl_matrix_free(fluxes);
    
}  // End of main routine



// **********************************************************
// ************  FUNCTION DEFINITIONS  **********************
// **********************************************************


// Spline interpolation of partition functions, if needed for this temperature.
// The current interpolated value will be stored in the isotope[] species
// objects and in the currentPF[] array for later use in constructing partition 
// function corrections for each reaction in the network that needs it.

void updatePF(){
    
    if(dopf && T9 > pfCut9){
        
        double pf_Iso[24];
        double pfNow;
        
        for (int i=0; i<ISOTOPES; i++){
            
            // Fill pf_Iso[] with pf table values for this isotope
            
            for(int k=0; k<PF; k++){
                pf_Iso[k] = isotope[i].getpf(k);
            }
            
            // Instantiate spline interpolator for partition function of this isotope
            
            SplineInterpolator interpolatePF = SplineInterpolator(PF, Tpf, pf_Iso);
            
            // Initialize spline interpolator for partition function of this isotope
            
            interpolatePF.spline((Tpf), pf_Iso, PF, PF);
            
            // Interpolate partition function for this isotope at present temperature
            
            pfNow = interpolatePF.splint(T9);
            
            // Store current value of pf in species objects and array currentPF[]
            
            isotope[i].setpfnow(pfNow);
            currentPF[i] = pfNow;

        }
    }
    
}  // End of updatePF()




// Set up headers for the plot files that will be written to during 
// the integration by toPlotNow(). Note that Lines beginning with # 
// are comments in gnuplot.

void plotFileSetup(){
    
    // Following arrays control which mass fractions are exported to plotting
    // file.  The entries in plotXlist[] are the species indices for the
    // isotopes in the network to be plotted. 
    
    for(int i=0; i<maxPlotIsotopes; i++){
        plotXlist[i] = i;
    }
    
    // Get length LX of array plotXlist holding the species indices for
    // isotopes that we will plot mass fraction X for.
    
    int LX = sizeof(plotXlist)/sizeof(plotXlist[0]);
    
    string str1 = "#    t     dt     |E|  |dE/dt| Asy  Equil sumX";
    string strflux = "\n#    t     dt   ";
    string app = "  ";
    string app1;
    string appflux;
    string Xstring = "  X(";
    string Fpstring = "F+(";
    string Fmstring = "F-(";
    string dFstring = "dF(";
    string iso;
    
    fprintf(pFileD,"\n\n");
    
    if(doASY){
        fprintf(plotfile1,"# ASY");
        fprintf(plotfile2,"# ASY");
        if(plotFluxes){fprintf(plotfile3,"# ASY");}
        fprintf(pFileD,"# ASY");
    } else {
        fprintf(plotfile1,"# QSS");
        fprintf(plotfile2,"# QSS");
        if(plotFluxes){fprintf(plotfile3,"# QSS");}
        fprintf(pFileD,"# QSS");
    }
    
    if(doPE){
        fprintf(plotfile1,"+PE");
        fprintf(plotfile2,"+PE");
        if(plotFluxes){fprintf(plotfile3,"+PE");}
        fprintf(pFileD,"+PE");
    } 
    
    if(dopf){
        fprintf(plotfile1," method (with partition functions). ");
        fprintf(plotfile2," method (with partition functions). ");
        if(plotFluxes){fprintf(plotfile3," method (with partition functions). ");}
        fprintf(pFileD,"+ method (with partition functions). ");
    } else {
        fprintf(plotfile1," method (no partition functions). ");
        fprintf(plotfile2," method (no partition functions). ");
        if(plotFluxes){fprintf(plotfile3," method (no partition functions). ");}
        fprintf(pFileD,"+ method (no partition functions). "); 
    }
    
    fprintf(plotfile1,"\n# All quantities except Asy, RG_PE, and sumX are \n");
    fprintf(plotfile1,"# log10(x). Log of absolute values for E and dE/dt, as\n");
    fprintf(plotfile1,"# they can be negative. UNITS: t and dt in s; E in erg; \n");
    fprintf(plotfile1,"# dE/dt in erg/g/s; all others dimensionless \n");
    fprintf(plotfile1,"#\n");
    
    // Write header for file pointed to by plotfile1
    
    for(int i=0; i<LX; i++){
        
        int indy = plotXlist[i];
        iso = to_string(indy);       // Use species index as label
        
        // Use instead isotope symbol as label
        
        // iso = charArrayToString(isoLabel[indy],isoLen);  // char array to string
        
        // The above conversion of a character array to a string leaves
        // 2 or 3 end of line characters after the isotope label in iso
        // that prevent printing correctly in fprint(pFile, stringToChar(str1))
        // below. It will print with ofstream as commented out below, but 
        // displaying 2-3 garbage symbols at the end. Remove those characters
        // with pop_back(). Unfortunately, since the isotope symbols are from 3 to 5 
        // characters long, 3 applications of pop_back removes the unwanted
        // characters but will in some cases remove the last character of the 
        // actual isotopic symbol. Thus switched to displaying the species number
        // rather than the isotopic symbol in the plot output.  More compact anyway.
        
        // iso.pop_back();  // Remove last character
        // iso.pop_back();
        // iso.pop_back();
        
        app.append(Xstring);
        app.append(iso);
        app.append(")    ");
        
    }
    
    str1.append(app);
    str1.append("\n");
    
    // Write header for plotfile1 -> plot1.data output. 
    
    fprintf(plotfile1,Utilities::stringToChar(str1));
    
    // Write header for plotfile2 -> plot2.data output
    
    string str2 = "#  t          dt   2/Rmin Reaction_Rmin       2/Rmax   Reaction_Rmax";
    str2 += ("     dt_FE    dt_EA  trial_dt   interpT   interpRho\n");
    fprintf(plotfile2,
            "\n# All double quantities are log10(x); rates in units of s^-1\n#\n");
    fprintf(plotfile2, Utilities::stringToChar(str2));
    
    // Write header for plotfile3 stream -> plot3.data (fluxes)
    
    if(plotFluxes){
        for(int i=0; i<LX; i++){
            int indy = plotXlist[i];
            iso = to_string(plotXlist[i]);  
            appflux.append(Fpstring);
            appflux.append(iso);
            appflux.append(")   ");
        }
        
        for(int i=0; i<LX; i++){
            int indy = plotXlist[i];
            iso = to_string(plotXlist[i]); 
            appflux.append(Fmstring);
            appflux.append(iso);
            appflux.append(")   ");
        }
        
        for(int i=0; i<LX; i++){
            iso = to_string(plotXlist[i]); 
            appflux.append(dFstring);
            appflux.append(iso);
            appflux.append(")   ");
        }
        
        strflux.append(appflux);
        fprintf(plotfile3,Utilities::stringToChar(strflux));
    }
    
    // Write header for plotfile4 stream -> plot4.data
    
    if(hydroProfile && plotHydroProfile)
        fprintf(plotfile4,"#  time          T          rho");
    
    // Write header for plotfile4 stream -> plot5.data
    
    fprintf(plotfile5,"# All quantities actual values (not log10) \n");
    fprintf(plotfile5,"#\n");
    
    fprintf(plotfile5, "# time(s)    dt(s)      dE(erg)  E(erg/g/s)    T(K)    rho(g/cm^3) 2/Rmax(s)");
    
    
    
}  // End of plotFileSetup()


// Function to write data to output files at each plot step.  Headers for 
// these data files are created in plotFileSetup() above. Data are output 
// at regular intervals of the integration to a file suitable for plotting. 
// Assumes the existence of a subdirectory gnu_out. May crash if this directory
// does not exist. Assuming gnuplot for plotting, but the output files
// are whitespace-delimited ascii, so any plotting program could be used
// to read it. 


void toPlotNow(){
    
    // Output to plotfile1 stream -> plot1.data. 
    
    if(EfromMasses){      //  EfromMasses should be true except for testing
        fprintf(plotfile1,"\n%+6.3f %+6.3f %6.3f %6.3f %5.3f %5.3f %5.3f",
            log10(t),log10(dt),log10( abs(ECON*EReleaseA)),
            log10( abs(ECON*dEReleaseA) ),(double)totalAsy/(double)ISOTOPES,
            (double)totalEquilRG/(double)numberRG,sumX);
    } else {
        fprintf(plotfile1,"\n%+6.3f %+6.3f %6.3f %6.3f %5.3f %5.3f %5.3f",
            log10(t),log10(dt),log10( abs(ECON*ERelease)),
            log10( abs(ECON*netdERelease) ),(double)totalAsy/(double)ISOTOPES,
            (double)totalEquilRG/(double)numberRG,sumX);
    }
    
    // Now add one data field for each X(i) in plotXlist[]. Add
    // GZ = 1e-24 to X in case it is identically zero since we are
    // taking log10.
    
    for(int j=0; j<maxPlotIsotopes; j++){
        fprintf(plotfile1," %5.3e", log10(X[j] + GZ));
    }
    
    // Output to plotfile2 stream -> plot2.data
    
    fprintf(plotfile2,
        "%7.4f %7.4f %7.4f %s %7.4f %s %7.4f %7.4f %7.4f %7.4e %7.4e\n",
            log10(t),log10(dt),log10(1.0/slowestCurrentRate),
            reacLabel[ slowestCurrentRateIndex],
            log10(2.0/fastestCurrentRate),
            reacLabel[fastestCurrentRateIndex],
            log10(dt_FE),log10(dt_EA),log10(dt),
            logTnow,logRhoNow
    );
    
    // Output to plotfile3 stream -> plot3.data (fluxes)
    
    if(plotFluxes){
        fprintf(plotfile3,"\n%+6.3f %+6.3f",log10(t),log10(dt));
        
        // Now add one data field for each FplusSumPlot. Add
        // GZ = 1e-24 to X in case it is identically zero since we are
        // taking the log.
        
        for(int j=0; j<maxPlotIsotopes; j++){
            fprintf(plotfile3," %5.3e",
                log10(abs( isotope[j].getfplus() + GZ) ));
        }
        
        for(int j=0; j<maxPlotIsotopes; j++){
            fprintf(plotfile3," %5.3e",
                log10(abs( isotope[j].getfminus() + GZ)));
        }
        
        for(int j=0; j<maxPlotIsotopes; j++){
            fprintf(plotfile3," %5.3e",
                    log10( abs(isotope[j].getfplus() - isotope[j].getfminus() 
                    + GZ) ));
        }
    }
        
    // Output to plotfile4 stream -> plot4.data.  These are the interpolated
    // temperature and density during the calculation as a function of time.
    // The output file hydroProfileInput.data contains the input hydro profile
    // from which the temperature and density were interpolated.
        
    if(hydroProfile && plotHydroProfile){
        fprintf(plotfile4,"\n%6.4e  %6.4e  %6.4e", log10(t),logTnow,logRhoNow);
    }
    
    // Output to plotfile5 stream -> plot5.data
    
    if(EfromMasses){      //  EfromMasses should be true except for testing
        fprintf(plotfile5,"\n%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e",
            t, dt, ECON*dEReleaseA, ECON*EReleaseA, T9*1e9, rho, 2.0/fastestCurrentRate);
    } else {
        fprintf(plotfile5,"\n%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e",
                t, dt, ECON*netdERelease, ECON*ERelease, T9*1e9, rho, 2.0/fastestCurrentRate);
    }
    
    // Flush the buffers holding output to plotting data files at each plot
    // output step so that plots can be made during a calculation if desired.
    
    fflush(plotfile1);
    fflush(plotfile2);
    fflush(plotfile3);
    fflush(plotfile4);
    fflush(plotfile5);
        
}  // End of toPlotNow()


// Function restoreBe8() to convert all 8-Be to alpha particles since lifetime 
// of 8-Be to decay to two alpha particles is short compared with typical 
// integration steps. Execute at end of integration step if both 8_Be and
// 4-He are present in the network.

void restoreBe8(){

    Y[indexAlpha] += 2.0*Y[indexBe8];
    X[indexAlpha] += 2.0*X[indexBe8];
    Y[indexBe8] = 0.0; 
    X[indexBe8] = 0.0; 
    
}


// Function showParameters() displays integration parameters.

void showParameters(){
    
    printf("\n\nIntegration using ");
    printf(Utilities::stringToChar(intMethod));
    
    if(dopf){
        printf(" (Partition function corrections applied for T9=%4.2f and above)", pfCut9);
    } else {
        printf(" (Partition function correction ignored)");
    }
    
    if(doPE)
    printf("\nPartial equilibrium constraints: useEquilY = %d useDevious = %d", 
        useEquilY, useDevious);
    
    string options = "Plot output for energy release ";
    if(EfromMasses){
        options += "dE(t) and E(t) computed from mass difference over dt";
    } else {
        options += "dE(t) and E(t) computed from Q-values and fluxes";
    }
    
    cout << "\n" << options;
   
    if(fixCNO){
        printf("\nCNO cycle fix used: startX_fixCNO = %5.2e", startX_fixCNO);
        if(totalTimeSteps > 0)
            printf(" Cycle fix started at t=%6.3e s", startX_fixCNO_time);
    }
    
    if(hydroProfile){
        printf("\nmassTol_asy=%5.2e massTol_PE=%5.2e\nsf=%5.2e equiTol=%5.2e equilTime=%5.2e",
               massTol_asy,massTol_asyPE,sf,equiTol,equilTime);
    } else {
        printf("\nT9=%6.3e (constant) rho=%6.3e (constant) massTol_asy=%5.2e massTol_PE=%5.2e\nsf=%5.2e equiTol=%5.2e equilTime=%5.2e",
               T9,rho,massTol_asy,massTol_asyPE,sf,equiTol,equilTime);
    }
    
    printf("\nmaxit=%d downbumper=%6.3f EpsA=%5.2e EpsR=%5.2e deviousMax=%5.3f",
           maxit, downbumper, EpsA, EpsR, deviousMax
    );
    cout << "\nNetwork: " << networkFile;
    cout << "  Rates: " << rateLibraryFile;
    if(hydroProfile) cout << "\nHydro profile: " << hydroFile;
    printf("\nIsotopes=%d Reactions=%d",ISOTOPES,SIZE);
    if(numberRG > 0) printf(" ReactionGroups=%d",numberRG);
    if(numberRG > 0) printf(" SingletRG=%d", numberSingletRG);
    if(totalTimeSteps > 0)
    printf("\nIntegration steps=%d totalIterations=%d IntegrationSteps_plotted=%d",
        totalTimeSteps,totalIterations,totalTimeSteps-totalTimeStepsZero);
    if(totalTimeSteps > 0) printf("\nMax dt iterations = %d at step %d (t=%5.3e)", 
        mostIterationsPerStep, maxIterationStep, maxIterationTime);
    if(totalTimeSteps > 0) Utilities::stopTimer();      // Stop timer and print integration time
    
    if(reacLibWarn){
        cout << "\n\n*******************************************************************";
        cout << "\n*** WARNING: Some temperatures were below library validity bound. *";
        printf("\n*** First occurrence for T=%5.3e K at t=%6.4e s.         *",
            warnT*1e9,warnTime);
        printf("\n*** Reaction rates were approximated by zero when T < 1e7 K.      *");
        cout << "\n*******************************************************************\n\n";
    }
    
    cout.flush();
}


// Function dE_halfstep() Computes and return dE/dt for a half-step
// in the integration based on Q-values.

double dE_halfstep(){
    
    double dE_half = 0.0;
    for(int i=0; i<SIZE; i++){
        dE_half += reaction[i].getQ() * reaction[i].getflux();
    }
    return dE_half;
    
}


/* Function restoreEquilibriumProg() to adjust populations at end of timestep when in PE
 * to correct for deviations from equilibrium during the timestep caused by reactions
 * not in equilibrium. This function sets the species abundances in each reaction group 
 * to their equilibrium value at the end of the numerical timestep. Then, the corrected 
 * abundances Y for all isotopes participating in partial equilibrium are averaged over
 * reaction groups if they participate in more than one reaction group. Finally, after 
 * Ys are updated by their average equilibrium values, all abundances are rescaled so 
 * that total nucleon number is conserved by the overall timestep (by requiring the sum
 * of the mass fractions X to equal one). Thus this function for restoring equilibrium 
 * does not require a matrix solution or Newton-Raphson iteration. It should scale 
 * approximately linearly with the network size for sparse networks. Contrast with the 
 * different more complicated option given in the restoreEquilibrium() method of the 
 * benchmark Java code. */


void restoreEquilibriumProg() {
    
    int countEquilRG = 0;
    int countEquilIsotopes = 0;
    sumX = Utilities::sumMassFractions();
    
    /* In general when we compute the equilibrium value of say alpha in the 
     * reaction group alpha+16O <-> 20Ne, we are computing it using non-equilibrium 
     * values of 16O and 20Ne (i.e.,their values will not be the values that they 
     * will have after this step). Add a while loop that permits iteration to try 
     * to fix this. Preliminary tests indicate it has small effect so set to one 
     * iteration for now. */
    
    int itcounter = 0;
    
    while (itcounter < 1) {
        
        countEquilRG = 0;          // Counts RG in equilibrim at this point
        countEquilIsotopes = 0;    // Counts isotopic species in equilibrated RGs
        
        itcounter ++;
        
        /* Compute equilibrium value of the Ys participating in equilibrium 
         * starting from the value of Y at the end of the numerical timestep,
         * presently stored in Y[i]. Do so by first setting Y0[i] to the current 
         * value of Y[i], which is the computed value at the END of the timestep. 
         * Then evolve that initial value to the corresponding equilibrium value 
         * algebraically by calculating the equilibrium value for that Y0[i] 
         * and setting Y[i] to it (a form of operator splitting within the 
         * network timestep). This work is done in evolveToEquilibrium(). */
        
        evolveToEquilibrium();
        
        // Inventory reaction groups in equilibrium counting the number in the
        // the variable countEquilRG, placing the index of each equilibrated 
        // RG in the array RGindy[], and counting the number of different isotopes
        // in at least one equilibrated RG in the variable countEquilIsotopes.

        int RGindy[numberRG] = {-1};  // Array to hold index of RG in equilibrium
        
        for (int i = 0; i < numberRG; i++) {
            
            // Process RG in equilibrium

            if ( RG[i].getisEquil() ) {
                
                // Store in RGindy[] the RG index for this equilibrated RG.
                
                RGindy[countEquilRG] = i;
                countEquilRG++;
                
                // Loop over the RG[i].niso isotopic species in this 
                // equilibrated RG,setting the species array entry
                // isotopeInEquil[] to true if a species is in a
                // RG in equilibrium.  Thus species entry in 
                // isotopeInEquil[] will be true if it is in at
                // least one equilibrated RG, and false if it is in
                // no equilibrated RGs.
                
                for (int j = 0; j < RG[i].getniso(); j++) {
                    
                    // Get species index for this isotope
                    
                    int speciesIndy = RG[i].getisoindex(j);
                    isotopeInEquil[speciesIndy] = true;
                    countEquilIsotopes++;

                }
            } 
        }
        
        // Loop over reaction groups in equilibrium and compute equilibrated
        // Y[] averaged over all reaction groups that are in equilibrium and
        // contain the isotope. (Generally each reaction is in only one 
        // reaction group but a given isotopic species can appear in many reaction 
        // groups.)
        
        int numberCases;    // Number of equilibrated RGs an isotope is in
        double Ysum;        // Sum of Ys for isotope in all equilbrated RGs
        int indy;           
        
        // Loop over all isotopes,checking for those in equilbrium in at 
        // least one RG
        
        for(int i=0; i<ISOTOPES; i++){
            
            numberCases = 0;
            Ysum = 0.0;
           
            // If isotope is in at least one equilibrated RG
            
            if (isotopeInEquil[i]) {
                
                // Find all equilibrated RGs that the isotope appears in
                
                for(int j=0; j<countEquilRG; j++){
                    indy = RGindy[j];           // index of RG in equil
                    if( RG[indy].getisEquil() ){
                        for(int k=0; k<RG[indy].getniso(); k++){
                            double addit = 0.0;
                            int RGisoindex = RG[indy].getisoindex(k);
                            if(i == RGisoindex) {
                                addit = RG[indy].getisoYeq(k);
                                Ysum += addit;
                                numberCases ++;
                            }

                        }
                    }
                }
            }
            
            // Store Y for each isotope averaged over all reaction groups in 
            // which it participates
            
            if(numberCases > 0) {
                Y[i] = Ysum/(double)numberCases;
                X[i] = Y[i]*(double)AA[i];
            }

        }

    }   // end while iteration loop
    
    
    // If we are computing Asy+PE, set up renormalization of all Ys 
    // and Xs so that this integration step conserves total particle 
    // number (sum of X = 1.0). Option to check mass fractions separately 
    // for isotopes participating in equilibrium and those not but don't 
    // presently use the fractions separately, only their sum. Whether
    // sum X normalized to one after PE evolution controlled by flag
    // normPEX (normally set to true).
    
    if(normPEX){
        
        sumXtot = Utilities::sumMassFractions();
        
        // Factor to enforce particle number conservation
        
        XcorrFac = 1.0 / sumXtot;
        
        // Loop over all isotopes and renormalize so sum X = 1
        // for this step.
        
        for(int i=0; i<ISOTOPES; i++){
            
            X[i] *= XcorrFac;
            Y[i] = X[i] / (double)AA[i];
            
        }
        
        // Total sumX now renormalized to one
        
        sumX = 1.0;
    
    }
    
}   // end restoreEquilibriumProg()


// ------------------------------------------------------------------------
// Function evolveToEquilibrium() to set abundances for reaction 
// groups in equilibrium to the current integrated value (end of timestep) 
// and then evolve the abundances algebraically to their equilibrium values.
// ------------------------------------------------------------------------

void evolveToEquilibrium() {
    
    // Set Y0 equal to current Y (at end of numerical step)
    // for all RG in equilibrium and recompute equilibrium

    for (int i = 0; i < numberRG; i++) {
        
        if ( RG[i].getisEquil() ) {
            
            // Loop over species in RG
            
            for (int j = 0; j < RG[i].getniso(); j++) {
                int indy = RG[i].getisoindex(j);
                Y0[indy] = Y[indy];
                RG[i].setisoY0(j,Y0[indy]);
            }
            
            // Compute equilibrium with new values of Y0
            
            RG[i].computeEquilibrium();
        }
    }
    
}


// ---------------------------------------------------------------
// Function isoIsInRG(int isoindex, rgindex) to return 
// true if isotope labeled by species index isoindex is in the RG 
// labeled by rgindex and false otherwise. Not presently used
// since ReactionGroup::speciesIsInRG(speciesIndex) duplicates
// this functionality. NOT PRESENTLY USED.
// --------------------------------------------------------------

bool isoIsInRG(int isoindex, int rgindex) {
    
    for (int j=0; j<RG[rgindex].getniso(); j++){
        if( RG[rgindex].getisoindex(j) == isoindex ) return true;
    }

    return false;
    
}


// Function getmaxdYdt() finds the maximum dY/dt for an 
// isotope in the network

void getmaxdYdt(){
    
    double ck;
    maxdYdt = 0.0;
    maxdYdtIndex = -1;
    for(int i=0; i<ISOTOPES; i++){
        ck = isotope[i].getdYdt();
        if( abs(ck) > abs(maxdYdt) ){
            maxdYdtIndex = i;
            maxdYdt = abs(ck);
        }
    }
    
}


// Function to read in hydro profile from file fileName if hydroProfile 
// is true. 

void readhydroProfile(char *fileName){
    
     char line[60];
     int numberEntries;
     int dummy;
     double Time;
     double Temp;
     double Rho;
    
     // Open a file for reading. Exit if the file can't be read
    
    FILE *fr; 
    fr = fopen (fileName,"r");
    if( fr == NULL ){
        printf ("*** File Input Error: No readable file named %s\n",fileName);
        exit(1);
    }
    
    int index = 0;
    int subIndex = -1;
    
    // Read lines until NULL encountered. Data lines can contain up to 60 characters.
    // The first line of the data file fileName contains labels and not data.
    // The second line contains the number of data points (lines) in the file, the
    // following lines contain the temperature and density data for each time.
    
    while(fgets(line,60,fr) != NULL){
        
        subIndex ++;
        
        if (subIndex == 1){  // Line containing number of entries
            
            sscanf(line,"%d",&numberEntries);
            
            if(numberEntries > maxHydroEntries){
                printf("\n\nERROR: Number of entries in hydro profile table of file %s (%d) ",
                    fileName,numberEntries);
                printf("\ninconsistent with spline arrays. "); 
                printf("Change the static constant maxHydroEntries ");
                printf("to a value \nof %d and recompile.\n\n",numberEntries);
                exit(1);
            }
            
        } else if (subIndex > 1) {    // Line containing t,T,rho data
            
            sscanf(line,"%d %lf %lf %lf",&dummy,&Time,&Temp,&Rho); 
            
            // Store hydro profile in three arrays. We will interpolate in log10
            // values, so store the time, temperature, and density from the
            // hydro values as log10 of their value.
            
            hydroTime[index] = log10(Time);   // Time
            hydroTemp[index] = log10(Temp);   // Temperature(time)
            hydroRho[index] = log10(Rho);     // Density(time)
            
            index ++;
        }
        
        hydroLines = index;            // # Data line entries read in
    
    }
    
    // Close the input file
    
    fclose(fr);
    
}    // End of function readhydroProfile(char *fileName)



/* Function readNetwork(char *filename) to read the network data file line 
 * by line, with the filename as argument. This file is expected to have 4 
 * lines per isotope with the line structure
 * 
 *	 isotopeSymbol A  Z  N  Y  MassExcess
 *	 pf00 pf01 pf02 pf03 pf04 pf05 pf06 pf07
 *	 pf10 pf11 pf12 pf13 pf14 pf15 pf16 pf17
 *	 pf20 pf21 pf22 pf23 pf24 pf25 pf26 pf27
 * 
 * where isotopeSymbol is an isotope label, A=Z+N is the atomic mass number,
 * Z is the proton number, N is the neutron number, Y is the initial abundance,
 * MassExcess is the mass excess in MeV, and the pf are 24 values of the partition 
 * function for that isotope at different values of the temperature that will form 
 * a table for interpolation in temperature. The assumed 24 values of the temperature 
 * for the partition function table are in array Tpf:
 * 
 * { 0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0,
 * 4.5,5.0,6.0,7.0,8.0,9.0,10.0 } 
 * 
 * in units of 10^9 K. All fields on a line are separated by a blank space and 
 * there is no whitespace in the isotopeSymbol. The type signature of these four 
 * lines corresponding to a single isotope is
 * 
 *	string int int int double double
 *	double double double double double double double double
 *	double double double double double double double double
 *	double double double double double double double double
 * 
 * Here is an example for two isotopes:
 * 
 * ca40 40 20 20 0.0 -34.846
 * 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 * 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 * 1.0 1.0 1.0 1.01 1.04 1.09 1.2 1.38
 * ti44 44 22 22 0.0 -37.548
 * 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 * 1.0 1.0 1.0 1.0 1.01 1.03 1.08 1.14
 * 1.23 1.35 1.49 1.85 2.35 3.01 3.86 4.94
 * 
 * A file with this format is written from the Java code to the file 
 * output/CUDAnetwork.inp using the Java stream toCUDAnet.
 *	
 */

void readNetwork (char *fileName) {
    
    char line[60];
    char isoSymbol[5];
    int z,n,a;
    double y,mass;
    double pf0,pf1,pf2,pf3,pf4,pf5,pf6,pf7;
    
    // Open a file for reading. Exit if the file can't be read
    
    FILE *fr; 
    fr = fopen (fileName,"r");
    if( fr == NULL ){
        printf ("*** File Input Error: No readable file named %s\n",fileName);
        exit(1) ;
    }
    
    // Read in the file line by line
    
    int isoIndex = -1;
    int isoSubIndex = 3;
    
    // Read lines until NULL encountered. Data lines can contain up to 60 characters.
    
    while(fgets(line, 60, fr) != NULL){
        
        isoSubIndex ++;
        
        if(isoSubIndex == 4){      // 1st line of data for first isotope
            
            isoSubIndex = 0;
            isoIndex ++;
            
            // Read 1st line
            
            sscanf (line,"%s %d %d %d %lf %lf",isoSymbol,&a,&z,&n,&y,&mass);
            
            // Write data in 1st line to the Species object isotope[]
            
            isotope[isoIndex].setisoIndex(isoIndex);
            isotope[isoIndex].setisoLabel(isoSymbol);
            isotope[isoIndex].setZ(z);
            isotope[isoIndex].setN(n);
            isotope[isoIndex].setA(a);
            isotope[isoIndex].setY(y);
            isotope[isoIndex].setM(mass);
            
        } else {             // line contains partition function entries
            
            // Scan and parse a partition function line. Each line contains 8
            // partition function values and there are three lines of partition
            // function values (so total of 24). The temperatures corresponding
            // to these 24 partition function values are given in the constant
            // array Tpf[PF] in the Species objects isotope[].
            
            sscanf (line,"%lf %lf %lf %lf %lf %lf %lf %lf", &pf0, &pf1, &pf2, &pf3,
                &pf4, &pf5, &pf6, &pf7);
            
            // Store the 24 partition function table values in Species object isotope[]
            
            int tin = isoSubIndex-1;
            
            // Set Species pointer to address of Species object isotope[isoIndex]
            
            SpeciesPtr = &isotope[isoIndex];
            
            SpeciesPtr->setpf(8*(tin), pf0);
            SpeciesPtr->setpf(8*(tin)+1, pf1);
            SpeciesPtr->setpf(8*(tin)+2, pf2);
            SpeciesPtr->setpf(8*(tin)+3, pf3);
            SpeciesPtr->setpf(8*(tin)+4, pf4);
            SpeciesPtr->setpf(8*(tin)+5, pf5);
            SpeciesPtr->setpf(8*(tin)+6, pf6);
            SpeciesPtr->setpf(8*(tin)+7, pf7);
        }
        
        // Normally numberSpecies = ISOTOPES,but numberSpecies counts the 
        // actual number of reactions read in.
        
        numberSpecies = isoIndex + 1;  

    }
    
//     Added static method to class Species to renorm all X to sum of one that is 
//     invoked after read-in of network. Now we start with sumX=1 to machine precision,
//     even if the conversion of the Ys read in to Xs does not give sumX = 1 to machine
//     precision.
    
    Species::renormAllX();
    
    // Close the file
    
    fclose(fr);
    
}    // End of function readNetwork (char *fileName)



/* Function readLibraryParams (char *fileName) to read rate parameter 
 * data file line by line, with fileName as argument. This file is expected 
 * to have one reaction per line with the line structure
 * 
 *    p0 p1 p2 p3 p4 p5 p6 reactionLabel
 * 
 * where the pn are the values of the 7 Reaclib parameters for a reaction,
 * reactionLabel is a label for the reaction that must contain no whitespace, and
 * all fields on a line are separated by a blank space.
 */

void readLibraryParams (char *fileName) {
    
    char line[120];
    char rlabel[LABELSIZE];
    double p0,p1,p2,p3,p4,p5,p6,q,sf;
    int i0,i1,i2,i3,i4,i5,i6,i7;
    int ii[6];
    double tempp[7];
    int tempZN[4] = {-1,-1,-1,-1};
    
    // Open a file for reading. Exit if the file can't be read.
    
    FILE *fr; 
    fr = fopen (fileName,"r");
    if( fr == NULL ){
        printf ("*** File Input Error: No readable file named %s\n",fileName);
        exit(1) ;
    }
    
    /* 
     * Read in the file line by line and parse into variables.  Each reaction corresponds
     * to 8 lines of input:
     * 
     *   reacLabel RGclass  RGmember reacClass #reac #prod  isEC isReverse Q prefac isForward
     *   p0        p1       p2       p3        p4    p5     p6   (reaction library constants)
     *   reac1Z    reac2Z   reac3Z                               (Z of reactants, #reac entries)
     *   reac1N    reac2N   reac3N                               (N of reactants, #reac entries) 
     *   prod1Z    prod2Z   prod3Z   prod4Z                      (Z of products, #prod entries)
     *   prod1N    prod2N   prod3N   prod4N                      (N of products, #prod entries)
     *   reac1Iso  reac2Iso reac3Iso                             (iso index reactants)
     *   prod1Iso  prod2Iso prod3Iso prod4Iso                    (iso index products)
     * 
     * where each entry is separated by a space, with NO WHITESPACE IN STRING reacLabel.
     */
    
    int n = -1;
    int subindex = -1;
    
    int newi0;
    
    // Read lines until NULL encountered. Lines can contain up to 120 characters.  In the
    // data file each reaction has 8 lines of entries.  The counter n holds the reaction number
    // (0 to SIZE) and the counter subindex holds the line number for the current reaction (0-7).  
    // The switch(subindex) determines which line we are reading (0-7) for a given reaction 
    // labeled by n.
    
    while(fgets(line,120,fr) != NULL) {
        
        subindex ++;
        
        switch(subindex) {
            
            case 0:   // 1st line
                
                n++;
                
                sscanf (line,"%s %d %d %d %d %d %d %d %lf %lf %d",
                    rlabel,&i0,&i1,&i2,&i3,&i4,&i5,&i6,&sf,&q,&i7);
                
                // Following converts any RGclass i0 = -1 produced by the Java code
                // to RG class 6 (which has only one reaction of the form 
                // a+b->c+d+e+f corresponding to ReacLib reaction class 7.)
                // Otherwise the -1 will be used as an index, causing a segfault.
                
                newi0 = i0;
                if(newi0 == -1) newi0 = 6;
                
                // Store in the Reaction class instance reaction[n]. First set a
                // pointer to the array of Reaction objects reaction[]
                
                ReactionPtr = &reaction[n];
                
                // Now point to various setter functions in the class Reaction and
                // use them to store the line of data just read in.
                
                ReactionPtr->setreacIndex(n);
                ReactionPtr->setreacString(rlabel);      // setter also fills reacLabel[][]
                ReactionPtr->setreacGroupClass(newi0);   // setter also fills RGclass[]
                ReactionPtr->setRGmemberIndex(i1);       // setter also fills RGMemberIndex[]
                ReactionPtr->setreacClass(i2);
                ReactionPtr->setnumberReactants(i3);     // setter also fills NumReactingSpecies[]
                ReactionPtr->setnumberProducts(i4);      // setter also fills NumProducts[]
                ReactionPtr->setisEC(i5);
                ReactionPtr->setisReverse(i6);
                ReactionPtr->setQ(q);
                ReactionPtr->setprefac(sf);
                ReactionPtr->setispeforward(i7);
                
                break;
                
            case 1:    // 2nd line (ReacLib rate coefficients)
                
                sscanf (line,"%lf %lf %lf %lf %lf %lf %lf",&p0,&p1,&p2,&p3,&p4,&p5,&p6);
                
                tempp[0] = p0;
                tempp[1] = p1;
                tempp[2] = p2;
                tempp[3] = p3;
                tempp[4] = p4;
                tempp[5] = p5;
                tempp[6] = p6;
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                ReactionPtr->setp(tempp);
                
                break;
                
            case 2:    // 3rd line (Z of reactants)
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantZ(tempZN);    // setter also changes reacZ[][]
                
                break;
                
            case 3:   // 4th line (N of reactants)
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantN(tempZN);   // setter also changes reacN[][]
                
                break;
                
            case 4:    // 5th line (Z of products)
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setproductZ(tempZN);    // setter also changes prodZ[][]
                
                break;
                
            case 5:    // 6th line (N of products)
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setproductN(tempZN);    // setter also changes prodN[][]
                
                break;
                
            case 6:    // 7th line (reactant species-vector index )
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantIndex(tempZN);   // setter also changes ReactantIndex[][]
                
                break;
                
            case 7:    // 8th line (product species-vector index)
                
                sscanf (line,"%d %d %d %d",&ii[0],&ii[1],&ii[2],&ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }

                ReactionPtr -> setproductIndex(tempZN);    // setter also changes ProductIndex[][]
                
                subindex = -1;
                
                break;
                
        }

    }
    
    // Normally numberReactions=SIZE, but numberReactions counts the actual number of
    // reactions read in.
    
    numberReactions = n+1;
    
    fprintf(pfnet,"\nSUMMARY: %d REACTIONS\n",numberReactions);
    
    // Close the file
    
    fclose(fr);           
    
}    // End of function readLibraryParams (char *fileName)



// Function writeNetwork() to send the network isotopes, mass excesses,
// and the entries in the partition function table for each isotope to the data file.

void writeNetwork() {

    fprintf(pfnet,"\n%d ISOTOPES IN NETWORK:\n\n",numberSpecies);
    fprintf(pfnet,"Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    
    for (int i=0; i<numberSpecies; i++){

        fprintf(pfnet ,"%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n", 
            i,isoLabel[i],AA[i],Z[i],N[i],
            Y[i],X[i],massExcess[i]);
        
    }
    
    // Write partition function table to output data file
    
    if(showDetails2) fprintf(pfnet,"\n\nPARTITION FUNCTION TABLE from Species object isotope[]:\n");
    
    if(showDetails2) fprintf(pfnet,"\nT9:  ");
    for(int k=0; k<24; k++){
        if(showDetails2) fprintf(pfnet,"%4.2f ", Tpf[k]);
    }
    
    if(showDetails2) fprintf(pfnet,"\nlgT: ");
    for(int k=0; k<24; k++){
        if(showDetails2) fprintf(pfnet,"%4.2f ", log10(Tpf[k]*1e9));
    }
    
    for(int i=0; i<ISOTOPES; i++){

        if(showDetails2) fprintf(pfnet,"\n");
        if(showDetails2) fprintf(pfnet,"%-5s", isotope[i].getLabel());
        
        for(int j=0; j<24; j++){ 
            if(showDetails2) fprintf(pfnet,"%4.2f ", isotope[i].getpf(j)); 
        }
        
    }

    if(showDetails2) fprintf(pfnet,"\n\n");
    
}   // End of function writeNetwork()



// Function writeRates(char *label) to print out all the rates. The 
// label can be used to distinguish cases if called more than once. NOT 
// PRESENTLY USED.

void  writeRates(char *label) {
    
    printf("\n\nCOMPUTED RATES (%s):\n\n",label);
    
    for (int i=0; i<numberReactions; i++){
        
        printf("%d %s rate=%6.3e Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Q=%5.3f Prefac=%6.3e Reactants=%d\n",
            i,reaction[i].getreacChar(), 
            Rate[i]/reaction[i].getprefac(),Rate[i],
            Y[reaction[i].getreactantIndex(0)],  
            Y[reaction[i].getreactantIndex(1)],
            Y[reaction[i].getreactantIndex(2)],
            reaction[i].getQ(),
            reaction[i].getprefac(),
            reaction[i].getnumberReactants());
        
    }
    
}    // End of function writeRates(char *label)


// Function setRG(int,int,int) that can be called from main or 
// any class to set some fields in the Reaction objects reaction[].

void setRG(int index,int RGclass,int RGindex) { 
    
    reaction[index].setrgindex(RGindex);
    
}


// Function that can be called from main or any class to create an array of 
// ReactionGroup objects RG[i] and to set some fields in the objects objects RG[].

void assignRG(){
    
    if(showDetails) fprintf(pFileD,"\n\nREACTIONS IN RGclass[]:\n");
    
    for(int m=0; m<SIZE; m++){
        if(showDetails) fprintf(pFileD,"\n%s RGclass[%d] = %d",reacLabel[m],m,RGclass[m]);
    }
    
    // Write out some fields for Reaction objects reaction[]
    
    if(showDetails) fprintf(pFileD,"\n\n\nSOME FIELDS FOR THE %d Reaction OBJECTS reaction[]:\n",SIZE);
    
    for(int i=0; i<SIZE; i++){
        
        if(showDetails) fprintf(pFileD,
            "\nreaction[%d]: %s RGclass=%d #reac=%d #prod=%d RGmemberIndex=%d RG=%d %s",
            i,reaction[i].getreacChar(),
            reaction[i].getreacGroupClass(),
            reaction[i].getnumberReactants(),
            reaction[i].getnumberProducts(),
            reaction[i].getRGmemberIndex(),
            reaction[i].getrgindex(),
            reacLabel[i]
        );
        
        int nummreac = reaction[i].getnumberReactants();
        int nummprod = reaction[i].getnumberProducts();
        
        // Write reactant symbols
        
        if(showDetails) fprintf(pFileD,"\nReaction=%d  REACTANTS: iso[0]=%s",
        i,isoLabel[reaction[i].getreactantIndex(0)]);
        
        if(nummreac > 1 && showDetails) fprintf(pFileD,"    iso[1]=%s",isoLabel[reaction[i].getreactantIndex(1)]);
        if(nummreac > 2 && showDetails) fprintf(pFileD," iso[2]=%s",isoLabel[reaction[i].getreactantIndex(2)]);
        
        // Write product Symbols
        
        if(showDetails) fprintf(pFileD,
        "  PRODUCTS: iso[%d]=%s",nummreac,isoLabel[reaction[i].getproductIndex(0)]);
        
        if(nummprod > 1 && showDetails) fprintf(pFileD," iso[%d]=%s",nummreac+1,isoLabel[reaction[i].getproductIndex(1)]);
        if(nummprod > 2 && showDetails) fprintf(pFileD," iso[%d]=%s",nummreac+2,isoLabel[reaction[i].getproductIndex(2)]);
        
        if(showDetails) fprintf(pFileD,"\n");
        
    }
    
    // Loop to create and populate ReactionGroup objects RG[]
    
    if(showDetails) fprintf(pFileD,"\n\nCREATING REACTION GROUPS RG[] AND POPULATING OBJECT FIELDS\n");
    
    for(int i=0; i<numberRG; i++){
        
        // Create array RG[i] of ReactionGroup objects. ReactionGroup(i) is the constructor
        // of the class ReactionGroup.
        
        RG[i] = ReactionGroup(i);
        RG[i].setnumberMemberReactions(RGnumberMembers[i]);
        
        // Set the reference reaction for the RG to be the first reaction in the RG
        // that has ifPEforward = true. This reference reaction will define the assumed
        // order of isotopes in the RG to be consistent with benchmark Java code.
        
        int reffer = RG[i].setrefreac();

        int rgindex = -1;
        int upper1;
        int upper2;
        int ck1;
        
        // Loop over reactions of the network,picking out the members of RG[i] by the
        // condition that i = RGindex[j], where RGindex[j] holds the RG index for a
        // given reaction.
        
        for(int j=0; j<SIZE; j++){   
            
            // Following if() condition picks from the list of reactions one that is
            // in the ReactionGroup object RG[i]
            
            if(RGindex[j] == i){
                
                rgindex ++;          // Index for member reactions in RG
                
                RG[i].setmemberReactions(rgindex,j);
                RG[i].setrgclass(RGclass[j]);
                RG[i].setisForward(rgindex,isPEforward[j] );
                
                ck1 = RG[i].getmemberReactions(rgindex);  //reacIndex of RG[i] member reaction
                
                RG[i].setnumberReactants(rgindex,reaction[ck1].getnumberReactants());
                RG[i].setnumberProducts(rgindex,reaction[ck1].getnumberProducts());
                RG[i].setreaclabel(rgindex,reaction[ck1].getreacString());
                
                upper1 = reaction[ck1].getnumberReactants();
                int rn = RG[i].getrefreac();
                int nrn = RG[i].getnumberReactants(rn);
                int ppp = RG[i].getmemberReactions(rn);
                
                // Loop over reactant isotopes within this reaction
                
                for(int k=0; k<nrn; k++){
                    int qqq = reaction[ppp].getreactantIndex(k);
                    RG[i].setisoindex(k,qqq);
                    RG[i].setreactantIsoIndex(k,qqq);
                    RG[i].setisoZ(k,reaction[ppp].getreactantZ(k));
                    RG[i].setisoN(k,reaction[ppp].getreactantN(k));
                    RG[i].setisoA(k,reaction[ppp].getreactantA(k));
                    RG[i].setisolabel(k,isoLabel[qqq]);
                }
                
                // Loop over product isotopes
                
                int nrn2 = RG[i].getnumberProducts(rn);
                upper2 = reaction[ck1].getnumberProducts();
                
                for(int k=0; k<nrn2; k++){
                    int qqq = reaction[ppp].getproductIndex(k);
                    RG[i].setisoindex(k+nrn,qqq);
                    RG[i].setproductIsoIndex(k,qqq);
                    RG[i].setisoZ(k+nrn,reaction[ppp].getproductZ(k));
                    RG[i].setisoN(k+nrn,reaction[ppp].getproductN(k));
                    RG[i].setisoA(k+nrn,reaction[ppp].getproductA(k));
                    RG[i].setisolabel(k+nrn,isoLabel[qqq]);
                }

                // Set the Ys in the RG
                
                int upk = RG[i].getnumberReactants(rn) + RG[i].getnumberProducts(rn);
                for(int k=0; k<upk; k++){
                    int yindex = RG[i].getisoindex(k);
                    RG[i].setisoY(k,Y[yindex]);
                }
            }
        }
        
        // Populate the species index array for this RG using reference reaction reffer
        
        int RGclassRef = RGclass[RG[i].getmemberReactions(reffer)];
        RG[i].setniso(RGclassRef);
        if(showDetails) fprintf(pFileD,"\nRG[%d]: refreac=%d RGclassRef=%d niso=%d Reactions=%d\n",
            i,RG[i].getrefreac(),RGclassRef,
            RG[i].getniso(),RG[i].getnumberMemberReactions() );
        
        // Set values of numberC for each RG
        
        switch ( RG[i].getrgclass() ){
            
            case 1:
                
                RG[i].setnumberC(1);
                break;
                
            case 2:
                
                RG[i].setnumberC(2);
                break;
                
            case 3:
                
                RG[i].setnumberC(3);
                break;
                
            case 4:
                
                RG[i].setnumberC(3);
                break;
                
            case 5:
                
                RG[i].setnumberC(4);
                break;
                
            // Bookkeeping. rgclass = 6 (F) can only contain one reaction from ReacLib class
            // 7 (a+b -> c + d + e + f), since there are no 4-body reactions in ReacLib
            // to serve as the inverse reaction.  But we will assume that all reactions
            // belong to a unique reaction group for bookkeeping, even if there is only one 
            // reaction in the reaction group so it can never equilibrate.
                
            case 6: 
                
                RG[i].setnumberC(0);
                break;
                
        }
    }
    
    // Summary of reaction groups
    
    for(int i=0; i<numberRG; i++){
        
        if(showDetails) fprintf(pFileD,"\n\nSummary: RG=%d", RG[i].getRGn());
        int numr = RG[i].getnumberMemberReactions();
        
        for(int j=0; j<numr; j++){
            
            int reacID = RG[i].getmemberReactions(j);
            if(showDetails) fprintf(pFileD,
                "\n%d %s iso[0]=%s iso[1]=%s iso[2]=%s iso[3]=%s",
                j,reacLabel[reacID],RG[i].getisolabel(0),
                RG[i].getisolabel(1),RG[i].getisolabel(2),
                RG[i].getisolabel(3)
            );
        }
    }
    
    // Check that this function has assigned isotope indices in each
    // reaction group consistent with the order used for partial equilibrium
    // in the original Java code.
    
    if(showDetails) fprintf(pFileD,"\n\n\nSUMMARY of order for isotopes for PE in ReactionGroup objects:");
    
    for(int i=0; i<numberRG; i++){
        
        int rn = RG[i].getrefreac();
        int upjj = RG[i].getnumberReactants(rn) + RG[i].getnumberProducts(rn);
        if(showDetails) fprintf(pFileD,
            "\n\nRG=%d  RGclass=%d %s Species Index:",
            i, RG[i].getrgclass(),
            reaction[RG[i].getmemberReactions(RG[i].getrefreac())].getreacGroupSymbol() 
        );
        
        for(int jj=0; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," iso[%d]=%d",jj,RG[i].getisoindex(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n     ");
        for(int jj=0; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," isolabel[%d]=%s",jj,RG[i].getisolabel(jj));
        }
        
        if(showDetails) fprintf(pFileD,
            "\n      REACTANTS: reactantIndex[0]=%d",RG[i].getreactantIsoIndex(0));
        for(int jj=1; jj<RG[i].getnumberReactants(rn); jj++){
            if(showDetails) fprintf(pFileD," reactantIndex[%d]=%d",jj,RG[i].getreactantIsoIndex(jj));
        }
        
        if(showDetails) fprintf(pFileD,
            "\n      PRODUCTS: productIndex[0]=%d",RG[i].getproductIsoIndex(0));
        for(int jj=1; jj<RG[i].getnumberProducts(rn); jj++){
            if(showDetails) fprintf(pFileD," productIndex[%d]=%d",jj,RG[i].getproductIsoIndex(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n      Z[0]=%d",RG[i].getisoZ(0));
        for(int jj=1; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," Z[%d]=%d",jj,RG[i].getisoZ(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n      N[0]=%d",RG[i].getisoN(0));
        for(int jj=1; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," N[%d]=%d",jj,RG[i].getisoN(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n      A[0]=%d",RG[i].getisoA(0));
        for(int jj=1; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," A[%d]=%d",jj,RG[i].getisoA(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n      isoY[0]=%8.5e",RG[i].getisoY(0));
        for(int jj=1; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," isoY[%d]=%8.5e",jj,RG[i].getisoY(jj));
        }
        
        if(showDetails) fprintf(pFileD,"\n      isoYeq[0]=%8.5e",RG[i].getisoYeq(0));
        for(int jj=1; jj<upjj; jj++){
            if(showDetails) fprintf(pFileD," isoYeq[%d]=%8.5e",jj,RG[i].getisoYeq(jj));
        }
    }
}       // End function assignRG()


// Function that can be called from main or any class to set field
// in the Species objects isotope[].

void setSpeciesfplus(int index,double fp){
    isotope[index].setfplus(fp);
}

// Function that can be called from main or any class to set field
// in the Species objects isotope[].

void setSpeciesfminus(int index,double fm){
    isotope[index].setfminus(fm);
    isotope[index].setkeff(fm/(isotope[index].getY0()+1.0e-30));
    keff[index] = isotope[index].getkeff();
}

// Function that can be called from main or any class to set field
// dydt in the Species objects isotope[].

void setSpeciesdYdt(int index,double dydt){
    isotope[index].setdYdt(dydt);
}


// Diagnostic function.  Not presently used.

void showY(){
    
    double sumXnow = 0.0;
    for(int i=0; i<ISOTOPES; i++){
        fprintf(pFileD,"\n   Y[%d]=%12.10e Y0[%d]=%12.10e X[%d]=%12.10e",
            i, Y[i], i, Y0[i], i, X[i] );
        sumXnow += X[i];
    }
        fprintf(pFileD,"\n   sumX=%12.10e sumXnow=%12.10e sumXlast=%12.10e",
        sumX, sumXnow, sumXlast
    );
}


// Function updateY0() that can be called from main or any class to
// update main Y0[] array,Species objects isotope[],and ReactionGroup
// objects RG[i] with current Y[] values.

void updateY0(){
    
    // Set Y0 in main array Y0[i] and in Species objects isotope[i]
    
    for(int i=0; i<ISOTOPES; i++){
        Y0[i] = Y[i];
        isotope[i].setY0(Y[i]);
    }
    
    // Set Y0 in ReactionGroup objects RG[i]
    
    for(int i=0; i<numberRG; i++){
        
        int jup = RG[i].getniso();
        
        if(jup < 0 || jup > 6){
            printf("\n*** ERROR: Invalid value niso=%d in updateY0() for RG=%d\n\n",jup,i);
            exit(-1);          // Exit program since something is corrupt
        }
        
        for(int j=0; j<jup; j++){
            int jj = RG[i].getisoindex(j);
            RG[i].setisoY0(j,Y[jj]);
        }
    }
    
}


// Function computeReactionFluxes() that can be called from main or another 
// class to set flux fields in the reaction[] objects.  Generally rates only
// need be recomputed if the temperature or density change, but fluxes must 
// be recomputed by calling this function any time either the rates or the 
// species populations change.

void computeReactionFluxes(){
    
    // Use functions of the Reaction class to compute fluxes.  We have instantiated
    // a set of Reaction objects in the array reaction[i],one entry for each
    // reaction in the network. Loop over this array and call the computeFlux()
    // function of Reaction on each object. Fluxes must be recomputed at each timestep
    // since they depend on the rates and the abundances. If temperature and density
    // are constant the rates won't change,but the fluxes generally will since
    // the abundances change even if the rates are constant as the network evolves.
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeFlux();
    }
    
    // Set the flux fields in the ReactionGroup objects RG[] from master flux
    // array Flux[].
    
    for(int i=0; i<numberRG; i++){
        RG[i].setRGfluxes();
    }
    
    // Call the static function Reaction::populateFplusFminus() to populate F+ and F-
    // for each isotope set up in setupFplusFminus() from master flux array computed
    // with setRGfluxes() above. We do it in this way because each flux is used in more than
    // one place but this way it is computed only once for each timestep.  This is efficient
    // because computing fluxes is the most time-consuming aspect of the algebraic explicit
    // algorithms.
    
    Reaction::populateFplusFminus();
    
    sumFplusFminus();     // Sum F+ and F- for each isotope
}


// Function sumFplusFminus() to sum the total F+ and 
// F- for each isotope. Moved here from original static function in class 
// Reaction to make diagnostics easier.

void sumFplusFminus(){
    
    int minny = 0;
    double accum;
    double dydt;
    
    // Loop over isotopes and reactions and sum F+ for each species
    
    for(int i=0; i < numberSpecies; i++){	
        
        // Sum F+ for each isotope
        
        if(i > 0) minny = FplusMax[i-1]+1;
        accum = 0.0;
        
        for(int j=minny; j<=FplusMax[i]; j++){
            accum += Fplus[j];
        }
        
        setSpeciesfplus(i,accum);      // Also sets FplusSum[i] = accum;
    }   
        
    // Loop over isotopes and reactions and sum F- for each species
    
    for(int i=0; i < numberSpecies; i++){
        
        minny = 0;
        if(i>0) minny = FminusMax[i-1]+1;
        accum = 0.0;
        
        for(int j=minny; j<=FminusMax[i]; j++){
            accum += Fminus[j];
        }
    
        setSpeciesfminus(i,accum);    // Also sets FminusSum[i] = accum and keff
        setSpeciesdYdt(i,FplusSum[i] - FminusSum[i]);
        
        dF[i] = FplusSum[i] - FminusSum[i];     // Store net flux
        
    }
    
}


// Function displayRGstatus() display to send to pfnet -> network.data the status 
// of each reaction group RG[].

void displayRGstatus(){
    
    fprintf(pfnet,"\n\n\nPARTIAL EQUILIBRIUM REACTION GROUPS (t=%6.4e, fracRGequil=%5.2f)", 
            t, eqFrac);
    
    string equilStatus;
    
    for (int i=0; i<numberRG; i++){
        if (RG[i].getisEquil()) {
            equilStatus = "TRUE ";
        } else {
            equilStatus = "FALSE";
        }
        fprintf(pfnet,
        "\n\nREACTION GROUP %d: isEquil=%s netflux=%6.3e thisDevious/maxDevious=%5.3f", 
        i, Utilities::stringToChar(equilStatus), RG[i].getnetflux(),
        RG[i].getthisDeviousRG()/deviousMax);
        
        double ftrue;
        double ftrueSum = 0;
        double rsign;
        
        for (int j=0; j<RG[i].getnumberMemberReactions(); j++){
            
            if(isPEforward[RG[i].getmemberReactions(j)]) {
                rsign = 1.0;
            } else {
                rsign = -1.0;
            }
            
            ftrue = reaction[RG[i].getmemberReactions(j)].getflux_true();
            ftrueSum += (rsign*ftrue);
            
            fprintf(pfnet,"\n    %s flux=%6.3e flux_true=%6.3e", 
                reacLabel[RG[i].getmemberReactions(j)], Flux[RG[i].getmemberReactions(j)],
                ftrue);
        }
        
        fprintf(pfnet,"\n    Sum flux_true=%6.3e", ftrueSum);
        
        fprintf(pfnet,"\nIsotope Equil Ratios:");
        for (int j=0; j<RG[i].getniso(); j++){
            fprintf(pfnet,"%4s(%5.3e) ", 
                isoLabel[RG[i].getisoindex(j)], RG[i].geteqRatio(j));
        }
        
        fprintf(pfnet,"\nIsotope Mass Fractions:");
        for (int j=0; j<RG[i].getniso(); j++){
            fprintf(pfnet,"%4s(%5.3e) ", 
                isoLabel[RG[i].getisoindex(j)], RG[i].getisoYeq(j)/RG[i].getisoA(j));
        }
        
    }
    
}


// ------------------------------------------------------------------------
// Method to return the current total mass (in MeV) of all isotopes in the
// network.
// ------------------------------------------------------------------------

double returnNetworkMass() {
    
    double sumM = 0;
    for (int i=0; i<ISOTOPES; i++){
        sumM += Y[i] * massExcess[i];
    }

    return sumM;
}


// --------------------------------------------------------------------------
// Method to calculate energy release by taking difference of network masses
// before and after timestep. Alternative to summing Q-values. Call at end
// of timestep. lastNetworkMass holds mass at beginning of timestep.
// --------------------------------------------------------------------------

void networkMassDifference() {
    
    networkMass = returnNetworkMass();
    netdEReleaseA = (lastNetworkMass - networkMass);
    dEReleaseA = netdEReleaseA / dt_saved;
    EReleaseA += netdEReleaseA;
    lastNetworkMass = networkMass;      // Initialize for next timestep
    
}

// Check to see whether network contains the main CNO cycle.  Returns
// true if all isotopes in the main CNO cycle are present, false otherwise.

bool checkForCNOCycle(){
    
    bool CNOisInNet = false;
    
    if( Utilities::isInNet(6, 6) &&    // 12C
        Utilities::isInNet(7, 6) &&    // 13N
        Utilities::isInNet(6, 7) &&    // 13C
        Utilities::isInNet(7, 7) &&    // 14N
        Utilities::isInNet(8, 7) &&    // 15O
        Utilities::isInNet(7, 8)       // 15N
    ) {
        
        CNOisInNet = true;
    }
    
    return CNOisInNet;
}

// Find the isotopic index for the CNO isotopes and the reaction index 
// for relevant reactions in the main CNO cycle

void indexCNOCycle(){
    
    // Find indices in the species array for the CNO isotopes
    
    index12C = Utilities::returnNetIndexZN(6, 6);
    index13N = Utilities::returnNetIndexZN(7, 6);
    index13C = Utilities::returnNetIndexZN(6, 7);
    index14N = Utilities::returnNetIndexZN(7, 7);
    index15O = Utilities::returnNetIndexZN(8, 7);
    index15N = Utilities::returnNetIndexZN(7, 8);
    index1H = Utilities::returnNetIndexZN(1, 0);
    index4He = Utilities::returnNetIndexZN(2, 2);
    
    printf("\n\nIndex CNOIsotopes: 12C=%d 13N=%d 13C=%d 14N=%d 15O=%d 15N=%d 1H=%d 4He=%d", 
        index12C, index13N, index13C, index14N, index15O, index15N, index1H, index4He);
    
    // Index relevant charged-particle reactions in main CNO cycle (each 
    // reaction rate has two components)

    bool resultCompare;
    
    // Find reaction indices for two components of p+n14-->o15
    
    char symbol1[] = "p+n14-->o15";
    
    for(int i=0; i<SIZE; i++){
        
        resultCompare = Utilities::compareTwoCharArrays(symbol1, reacLabel[i]);
        
        if(resultCompare){
            
            if(index14N_pgamma[0] == -1){
                index14N_pgamma[0] = i; 
            } else {
                index14N_pgamma[1] = i;
            }
        }
        
    }
    printf("\nReaction: %s index[0]=%d index[1]=%d", 
           symbol1, index14N_pgamma[0], index14N_pgamma[1]);
    
    // Find reaction indices for two components of p+c13-->n14
    
    char symbol2[] = "p+c13-->n14";
    
    for(int i=0; i<SIZE; i++){
        
        resultCompare = Utilities::compareTwoCharArrays(symbol2, reacLabel[i]);
        
        if(resultCompare){
            
            if(index13C_pgamma[0] == -1){
                index13C_pgamma[0] = i; 
            } else {
                index13C_pgamma[1] = i;
            }
        }
        
    }
    printf("\nReaction: %s index[0]=%d index[1]=%d", 
           symbol2, index13C_pgamma[0], index13C_pgamma[1]);
    
    // Find reaction indices for two components of p+n15-->he4+c12
    
    char symbol3[] = "p+n15-->he4+c12";
    
    for(int i=0; i<SIZE; i++){
        
        resultCompare = Utilities::compareTwoCharArrays(symbol3, reacLabel[i]);
        
        if(resultCompare){
            
            if(index15N_pgamma[0] == -1){
                index15N_pgamma[0] = i; 
            } else {
                index15N_pgamma[1] = i;
            }
        }
        
    }
    printf("\nReaction: %s index[0]=%d index[1]=%d", 
           symbol3, index15N_pgamma[0], index15N_pgamma[1]);
    
}


// Try removing stiffness associated with beta decay in main CNO cycle by 
// computing the 15N and 13C populations at equilibrium from the 14N population,
// by requiring that in the closed cycle the rates must be equal around the 
// cycle at equilibrium. The boolean fixCNO controls whether this fix
// is applied if X[H] < startX_fixCNO.

void correctCNOCycle(){
    
    // 15N population
    
    double Y15Nfac = (Rate[index14N_pgamma[0]] + Rate[index14N_pgamma[1]])
        /(Rate[index15N_pgamma[0]] + Rate[index15N_pgamma[1]]);
        
    double Y15Ncorr = Y15Nfac * Y[index14N];   // Corrected abundance of 15N
    
    // 13C population
    
    double Y13Cfac = (Rate[index14N_pgamma[0]] + Rate[index14N_pgamma[1]])
    /(Rate[index13C_pgamma[0]] + Rate[index13C_pgamma[1]]);
    
    double Y13Ccorr = Y13Cfac * Y[index14N];   // Corrected abundance of 13C

    // Keep track of when CNO corrections started
    
    if(!fixingCNO_now){
        startX_fixCNO_time = t;
        fixingCNO_now = true;
        printf("\n***** Start CNO fix: time=%5.3e", t);
    }
    
    // Correct 15N abundances and mass fractions
    
    Y[index15N] = Y15Ncorr;
    X[index15N] = Y[index15N] * 15;
    
    // Correct 13C abundances and mass fractions
    
    Y[index13C] = Y13Ccorr;
    X[index13C] = Y[index13C] * 13;
    
}


