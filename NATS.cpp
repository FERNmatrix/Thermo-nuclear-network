/*
 * Code to implement explicit algebraic integration of astrophysical thermonuclear networks.
 * Execution assuming use of Fedora Linux and GCC compiler: Compile with
 * 
 *     gcc NATS.cpp -o NATS -lgsl -lgslcblas -lm -lstdc++
 * 
 * Resulting compiled code can be executed with
 * 
 *     ./explicitMatrix  | tee temp.txt
 * 
 * where | tee temp.txt is unix shell script outputting to screen and also piped to a file temp.txt. 
 * Execution for other Linux systems, or Mac or PC, will depend on the C/C++ compiler installed on 
 * your machine but should be similar.
 *  
 * 
 * To set up specific calculation:
 * 
 * 1. Change ISOTOPES and SIZE
 * 2. Change two input files for networkFile and rateLibraryFile
 * 3. Change doAsy, doQSS, and doPE to choose Asy, Asy+PE, QSS, QSS+PE options
 * 4. Change control parameters like stop_time, massTol, ...
 * 5. Change plot output mask plotXlist[]
 * 6. Change values of T9_start and rho_start
 *
 * 
 * AUTHORS:
 * ---------------
 * Nick Brey
 * Ashton DeRousse
 * Adam Cole
 * Mike Guidry
 * ----------------
 */


#include <stdio.h>
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

using namespace std;
using std::string;

// NOTE:  For different networks you must change the values of array sizes SIZE and ISOTOPES, and input
//        filenames rateLibraryFile[] and networkFile[] to appropriate values. For now these are hardwired
//        definitions but eventually we probably want to read those quantities in as data.

// SIZE defines the number of reactions to be calculated.  Need min of SIZE=4395 for 365-element network, 
// 1603 for 150-isotope network, 2234 for the 194-isotope network, 3176 for the 268-isotope network,
// 1566 for the nova134 network, 48 for the alpha network, 8 for 3-alpha 4He-12C-16O network, and 28 for 
// the 7-isotope pp network. ISOTOPES defines the number of isotopes in each network (for example, ISOTOPES=7
// for the 7-isotope pp-chain network.  These sizes are hardwired for now but eventually we may want to read
// them in and assign them dynamically.

#define ISOTOPES 16                 // Max isotopes in network (e.g. 16 for alpha network, 7 in pp, 3 in triple alpha)
#define SIZE 48                        // Max number of reactions (e.g. 48 for alpha network, 28 in pp, 8 in triple alpha)
#define plotSteps 300                 // Number of plot output steps

#define LABELSIZE 35                  // Max size of reaction string a+b>c in characters
#define PF 24                         // Number entries partition function table for isotopes
#define THIRD 0.333333333333333
#define TWOTHIRD 0.66666666666667
#define ECON 9.5768e17                // Convert MeV/nucleon/s to erg/g/s
#define LOG10 0.434294481903251       // Conversion natural log to log10
#define MEV 931.494                   // Conversion of amu to MeV

#define unitd static_cast<double>(1.0)  // Constant double equal to 1
#define zerod static_cast<double>(0.0)  // Constant double equal to 0

// Define some CPU timing utilities. Usage:
//     START_CPU;
//     ... code to be timed ...
//     STOP_CPU;
//     PRINT_CPU
// These may be used to time CPU processes.

clock_t startCPU, stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("Timer: %7.4e ms used", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU (fprintf(pFile, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU2 (fprintf(pFile2, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPUD (fprintf(pFileD, "in %g seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define PRINT_CPU_TEST (printf("\nTimer Test: %g ms used by CPU\n", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));

// File pointer for data read-in

FILE *fr;            

// Filename for input rates library data. The file rateLibrary.data output by the Java code through 
// the stream toRateData has the expected format for this file.  Standard test cases: 
// rateLibrary_alpha.data, rateLibrary_150.data, rateLibrary_365.data, rateLibrary_nova134.data,
// rateLibrary_3alpha.data, rateLibrary_pp.data.

char rateLibraryFile[] = "data/rateLibrary_alpha.data";

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format for 
// this file. Standard test cases: CUDAnet_alphasolar.inp, CUDAnet_150solar.inp,
// CUDAnet_365solar.inp, CUDAnet_nova134.inp, CUDAnet_3alpha.inp, CUDAnet_pp.inp.

char networkFile[] = "data/CUDAnet_alpha.inp";

// File pointer for diagnostics output

FILE *pFileD;

// Control diagnostic printout of details (true=1 to print, false=0 to suppress)

static const int displayInput = true;
static const int showParsing = false;
static const int showFparsing = false;
static const int showFluxCalc = false;
static const int showRVdetails = false;
static const int showRGsorting = false;
static const int showAsyTest = false;
static const int showPlotSteps = false;
static const bool showAddRemove = false; 
static const bool showRestoreEq = false;
static const bool plotFluxes = false;
static const bool diagnose1 = false;
static const bool diagnose2 = false;
static const bool diagnose_dt = false;
static const bool diagnoseQSS = false;


// Function signatures in main:

void devcheck(int);
void readLibraryParams(char *);
void readNetwork(char *);
void writeNetwork(void);
void writeRates(char *);
void writeAbundances(void);
void setRG(int, int, int);
void setSpeciesfplus(int, double);
void setSpeciesfminus(int, double);
void setSpeciesdYdt(int, double);
void assignRG(void);
void plotOutput(void);
void getmaxdYdt(void);
void restoreEquilibriumProg(void);
void evolveToEquilibrium(void);
bool isoIsInRG(int, int);
void setReactionFluxes();

// Control which explicit algebraic approximations are used. Eventually
// this should be set from a data file. To use asymptotic set doASY true
// (which toggles doQSS to false). To use quasi-steady-state (QSS), set 
// doASY false (which toggles doQSS to true). doPE can be true or false 
// with either Asymptotic or QSS.

bool doASY = true;           // Whether to use asymptotic approximation
bool doQSS = !doASY;          // Whether to use QSS approximation 
bool doPE = true;             // Implement partial equilibrium also

// Temperature and density variables. Temperature and density can be
// either constant, or read from a hydro profile as a function of time.

double T9_start;              // Start temperature in units of 10^9 K
double rho_start;             // Start density in units of g/cm^3
double T9;                    // Current temperature in units of 10^9 K
double rho;                   // Current density in units of g/cm^3
bool constant_T9 = true;      // Whether temperature constant in integration
bool constant_rho = true;     // Whether density constant in integration

// Energy variables (from Q values)

double ERelease;              // Total energy released
double dERelease;             // Energy released per unit time

// Partition function controls. If dopf = true, reaction rates are
// corrected by temperature-dependent partition functions.  However
// partition function factors differ from 1 only at high temperature
// so we only implement partition function corrections if T9 > pfCut9,
// where pfCut9 is a cutoff temperature in units of T9. Typically we 
// will choose dopf = true and pfCut9 = 1.0.

bool dopf = true;
double pfCut9 = 1.0;

// Array to hold whether given species satisfies asymptotic condition
// True (1) if asyptotic; else false (0).

bool isAsy[ISOTOPES];
double asycheck;              // Species asymptotic if asycheck > 1.0

// Whether isotope part of any RG in partial equilibrium this timestep

bool isotopeInEquil [ISOTOPES]; 

// isotopeInEquil[] from last timestep

bool isotopeInEquilLast [ISOTOPES]; 

// Force a constant timestep constant_dt for testing purposes by
// setting constantTimestep=true.  Normally constantTimestep=false
// for adaptive timestepping.

bool constantTimestep = false;    // Adaptible timestep if false
double constant_dt = 1.1e-9;      // Value of constant timestep

// Integration time data.  The variables start_time and stop_time 
// define the range of integration (all time units in seconds),
// and dt_start sets the initial integration timestep. In an operator-split 
// coupling  start_time will be ~0, stop_time will correspond to the length
// of the hydro timestep and dt_init will likely be something like the 
// last timestep of the previous network integration (for the preceding 
// hydro timestep). Here we hardwire them for testing purposes.
// The variable startplot_time allows the plotting interval output
// in gnu_out/gnufile.data to be a subset of the full integration interval. 
// Generally, startplot_time > start_time.  By default the stop time for
// plotting is the same as the stop time for integration, stop_time.

double start_time = 1.0e-20;           // Start time for integration
double logStart = log10(start_time);   // Base 10 log start time
double startplot_time = 1.0e-11;       // Start time for plot output
double stop_time = 1.0e-2;             // Stop time for integration
double logStop = log10(stop_time);     // Base-10 log stop time
double dt_start = 0.01*start_time;     // Initial value of integration dt

double massTol = 1e-7;               // Timestep tolerance parameter (1.0e-7)
double SF = 7.3e-4;                    // Timestep agressiveness factor (7.3e-4)

// Time to begin trying to impose partial equilibrium if doPE=true. Hardwired but 
// eventually should be determined by the program.  In the Java version this was sometimes
// needed because starting PE test too early could lead to bad results.  This is 
// probably an error in the Java version, since if operating properly nothing should
// be changed at a timestep if nothing satisfies PE condition.  Thus, we should not need
// this in a final version for stability, but it might still be useful since early in
// a calculation typically nothing satisfies PE, so checking for it is a waste of time.
// On the other hand, check should not be costly.

double equilibrateTime = 1.0e-6;   // Begin checking for PE
double equiTol = 0.01;             // Tolerance for checking whether Ys in RG in equil

double deviousMax = 0.5;      // Max allowed deviation from equil k ratio in timestep
double deviousMin = 0.1;      // Min allowed deviation from equil k ratio in timestep
double mostDevious;           // Largest current deviation of k ratio from equil
int mostDeviousIndex;         // Index of RG with mostDevious


// Threshold abundance for imposing equil in reactions.  There may be numerical
// issues if the PE algorithm is imposed for very small abundances early in
// the calculation.

double Ythresh = 0.0;

double dt;                           // Current integration timestep
double dtgrow = 1.02;                // growth factor of dt
double dtdec = 0.98;                 // decrementing factor of dt
// double const uptol = 1e-7;        // upper tolerance condition
// double const lowtol = 1e-10;      // lower tolerance condition
// double tolC = 1e-7;               // %difference tolerance
double t;                            // Current time in integration
int totalTimeSteps;                  // Number of integration timesteps taken
double deltaTime;                    // dt for current integration step
int totalAsy;                        // Total number of asymptotic isotopes

double dtLast;                       // Last timestep
bool isValidUpdate;                  // Whether timestep accepted

double sumX;                         // Sum of mass fractions X(i).  Should be 1.0.
double sumXlast;                     // sumX from last timestep
double diffX;                        // sumX - 1.0

double maxdYdt;                      // Maximum current dY/dt in network
int maxdYdtIndex;                    // Isotopic index of species with max dY/dt

double Rate[SIZE];                   // Reaction rate from Reactions::computeRate()
double Flux[SIZE];                   // Flux from Reactions::computeFlux()


// --- Species data in following arrays also contained in fields of class Species

int Z[ISOTOPES];                 // Array holding Z values for isotopes
int N[ISOTOPES];                 // Array holding N values for isotopes
int AA[ISOTOPES];                // Array holding A values for isotopes
double Y[ISOTOPES];              // Array holding current abundances Y for isotopes
double Y0[ISOTOPES];             // Array holding abundances at beginning of timestep
double X[ISOTOPES];              // Array holding mass fractions X for isotopes
double massExcess[ISOTOPES];     // Array holding mass excesses for isotopes
char isoLabel[ISOTOPES][5];      // Isotope labels (max 5 characters; e.g. 238pu)
double dYDt[ISOTOPES];           // Rate of change for Y

// -----------


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

// -----------


int numberSpecies;                // Actual # species in network (usually = ISOTOPES)
int numberReactions;              // Actual # reactions in network (usually = SIZE)

// Array with entries +1 if a reaction increases the population of the isotope 
// (contributes to F+), -1 if it decreases it (contributes to F-) and 0 if the 
// reaction does not change the population of the isotope. This array is populated 
// in the function ReactionVector::parseF().  It is characteristic of the structure 
// of the network and thus has to be calculated only once for a given network.

int reacMask[ISOTOPES][SIZE]; 

// Define an array rv[] and corresponding pointers that will hold GSL vectors 
// corresponding to the reaction vectors for the system.  This will be implemented 
// in the function makeReactionVectors() of the class ReactionVector.

gsl_vector rv[SIZE];   // Array of type gsl_vector to hold GSL vectors
gsl_vector *rvPt;      // Pointer to rv[] array

// Total number of F+ and F- terms in the network

int totalFplus = 0;
int totalFminus = 0;

// Arrays to hold non-zero fluxes in the network. Corresponding memory will be allocated dynamically 
// below with malloc

double* Fplus;           // Dynamical 1D array for non-zero F+ (Dim totalFplus)
double* Fminus;          // Dynamical 1D array for non-zero F- (Dim totalFminus)

double FplusZero[ISOTOPES];  // Used in QSS approximation

double keff[ISOTOPES];        // Effective decay constant for given isotope
double keffZero[ISOTOPES];    // Used in QSS approximation

// Arrays to hold number species factors for F+ and F- arrays. For example, for 12C+12C ->
// 4He + 20Ne the species factor is 1 for 4He and 20Ne diff. equation terms but 2 for
// 12C diff equation term (since the reaction involves one 4He and one 20Ne, but two 12C.

double* FplusFac;     // Dynamically allocated 1D array for number species factor in F+ terms
double* FminusFac;    // Dynamically allocated 1D array for number species factor in F- terms

double* FplusSum;     // Sum of F+ for each isotope
double* FminusSum;    // Sum of F- for each isotope

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
 * master flux array Flux[SIZE]. For example, MapFplus[0] will hold the value of
 * the index i in Flux[i] that corresponds to the reaction that generates Fplus[0].
 * These will be allocated dynamically below.  In general a given reaction flux appears
 * in more than one Fplus or Fminus array.  So we compute the flux for each reaction only
 * once for a given temperature and density, and then use this mapping to assign a given
 * flux to its occurrence in multiple Fplus and Fminus entries.
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
int RGnumberMembers[SIZE];    // # members each RG; set in class ReactionVectors

// Define array to hold the reaction group index for each reaction. There are n reaction
// groups in the network and each reaction belongs to one reaction group.  RGindex[m] is
// the index (0, 1, ... n) of the reaction group to which reaction m belongs. 
// This array is populated by the function ReactionVector::sortReactionGroups()

int RGindex[SIZE];

// Could save a little memory by assigning dynamically to size [numberRG][ISOTOPES] 
// once numberRG determined

bool RGisoMembers[SIZE][ISOTOPES];

double Yminner;             // Current minimum Y in reaction group
double mineqcheck;          // Current minimum value of eqcheck in reaction group
double maxeqcheck;          // Current max value of eqcheck in reaction group

bool reacIsActive[SIZE];    // False if reaction has been removed by PE

int totalEquilReactions;    // Total equilibrated reactions for isotope
int totalEquilRG;           // Total equilibrated reaction groups

gsl_matrix *fluxes;
gsl_vector *abundances;

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


// -------------------------------------------------------------------------
// Arrays to hold output quantities at plot timesteps. The variable
// plotSteps is hardwired above, but eventually should be input
// at time of execution.
// -------------------------------------------------------------------------

// Target plot times for plot steps. Computed from value of plotSteps in
// Utilities::log10Spacing().

double plotTimeTargets[plotSteps];

double tplot[plotSteps];                     // Actual time for plot step
double dtplot[plotSteps];                    // dt for plot step
double Xplot[ISOTOPES][plotSteps];           // Mass fractions X
double sumXplot[plotSteps];                  // Sum of mass fractions
//double diffXplot[plotSteps];               // Deviation of sumX from 1.0
int numAsyplot[plotSteps];                   // Number asymptotic species
int numRG_PEplot[plotSteps];                 // Number RG in PE
double EReleasePlot[plotSteps];              // Integrated energy release
double dEReleasePlot[plotSteps];             // Differential energy release
double fastestRatePlot[plotSteps];           // Fastest reaction rate at time t
int fastestRateIndexPlot[plotSteps];         // Reaction index fastest rate
double slowestRatePlot[plotSteps];           // Slowest reaction rate at time t
int slowestRateIndexPlot[plotSteps];         // Reaction index slowest rate
double FplusSumPlot[ISOTOPES][plotSteps];    // FplusSum
double FminusSumPlot[ISOTOPES][plotSteps];   // FplusSum

// Following control which mass fractions are exported to the plotting
// file.  The entries in plotXlist[] are the species indices for the
// isotopes in the network to be plotted.

//int plotXlist[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};   //alpha
int plotXlist[] = {1, 2, 3};    // 3-alpha
int LX;                         // Length of plotXlist array




//----------------CLASS DEFINITIONS ----------------


/* Class Utilities to hold utility useful utility functions.  Functions are
 * declared static so that they can be invoked without having to instantiate
 * objects of type Utilities.  For example, Utilities::returnNetIndexZN (Z, N).
 * We will often have other classes inherit from Utilities so that they can
 * access its static functions.  This class definition must precede definitions of 
 * classes that inherit from it.
 */


class Utilities{
    
    private:
    
    public:
        
        // Static function Utilities::showTime() to return date and local time as a 
        // character array
        
        static char* showTime(){
            time_t now = time(0);         // Current date/time
            char* tnow = ctime(&now);     // convert to string
            return tnow;
        }
        
        // -------------------------------------------------------------------------
        // Static function Utilities::interpolate_T(t) to find an interpolated T9
        // as a function of time if constant_T9 = false.
        // -------------------------------------------------------------------------
        
        static double interpolate_T(double t){
            
            // Will call spline interpolator in hydro profile table to return 
            // T9 at this value of time t.  For now, we just return T9_start.
            
            double temp;            // Interpolated temperature
            temp = T9_start;        // Temporary for testing
            
            printf("\n\n**** Temperature interpolation T9(%7.4e) = %6.4f ****",
                t, temp
            );
            
            return temp;
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::interpolate_rho(t) to find an interpolated 
        // density rho as a function of time if constant_rho = false.
        // -------------------------------------------------------------------------
        
        static double interpolate_rho(double t){
            
            // Will call spline interpolator in hydro profile table to return 
            // rho at this value of time t.  For now, we just return rho_start.
            
            double rhonow;             // Interpolated density
            rhonow = rho_start;        // Temporary for testing
            
            printf("\n**** Density interpolation rho(%7.4e) = %7.4e ****\n",
                   t, rhonow
            );
            
            return rhonow;
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::log10Spacing() to find num equal log10 
        // spacings between two numbers start and stop. The equally log-spaced
        // numbers are placed in the array v passed with the pointer *v.
        // -------------------------------------------------------------------------
        
        static void log10Spacing(double start, double stop, int num, double *v){
            
            double logtmin = log10(start);
            double logtmax = log10(stop);
            double tempsum = logtmin;
            double expofac = (logtmax - logtmin) / (double) num;
            v[0] = pow(10, logtmin);
            for(int i=0; i<num; i++){
                tempsum += expofac;
                v[i] = pow(10, tempsum);
            }
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::plotOutput() to output data at regular 
        // intervals of the integration to a file suitable for plotting. Assumes
        // the existence of a subdirectory gnu_out. Will crash if this directory
        // does not exist. Assuming gnuplot for plotting, but the output file
        // is whitespace-delimited ascii, so any plotting program could be used
        // to read it. Lines beginning with # are comments in gnuplot.
        // -------------------------------------------------------------------------
        
        static void plotOutput(){
            
            // Open files for ascii output. Assumes that the subdirectory
            // gnu_out already exists. If it doesn't, will compile but
            // crash when executed.
            
            FILE * pFile;
            pFile = fopen("gnu_out/gnufile.data","w");
            
            FILE * pFile2;
            pFile2 = fopen("gnu_out/gnufile2.data","w");
            
            FILE * pFile3;
            pFile3 = fopen("gnu_out/gnufileFlux.data","w");
            
            // Get length of array plotXlist holding the species indices for isotopes
            // that we will plot mass fraction X for.
            
            LX = sizeof(plotXlist)/sizeof(plotXlist[0]);
            
            string str1 = "#    t     dt     |E|  |dE/dt| Asy  Equil  sumX"; //add diffX to keep track of it at every step
            string strflux = "\n#    t     dt   ";
            string app = "  ";
            string app1;
            string appflux;
            string Xstring = "X(";
            string Fpstring = "F+(";
            string Fmstring = "F-(";
            string iso;
            
            fprintf(pFileD, "\n\n");
            
            if(doASY){
                fprintf(pFile, "# ASY");
                fprintf(pFile2, "# ASY");
                if(plotFluxes){fprintf(pFile3, "# ASY");}
                fprintf(pFileD, "# ASY");
            } else {
                fprintf(pFile, "# QSS");
                fprintf(pFile2, "# QSS");
                if(plotFluxes){fprintf(pFile3, "# QSS");}
                fprintf(pFileD, "# QSS");
            }
            
            if(doPE){
                fprintf(pFile, "+PE");
                fprintf(pFile2, "+PE");
                if(plotFluxes){fprintf(pFile3, "+PE");}
                fprintf(pFileD, "+PE");
            } 
            
            if(dopf){
                fprintf(pFile, " method (with partition functions): ");
                fprintf(pFile2, " method (with partition functions): ");
                if(plotFluxes){fprintf(pFile3, " method (with partition functions): ");}
                fprintf(pFileD, "+ method (with partition functions): ");
            } else {
                fprintf(pFile, " method (no partition functions): ");
                fprintf(pFile2, " method (no partition functions): ");
                if(plotFluxes){fprintf(pFile3, " method (no partition functions): ");}
                fprintf(pFileD, "+ method (no partition functions): "); 
            }
            
            fprintf(pFile, "%d integration steps ", totalTimeSteps);
            fprintf(pFile2, "%d integration steps ", totalTimeSteps);
            if(plotFluxes) fprintf(pFile3, "%d integration steps ", totalTimeSteps);
            fprintf(pFileD, "%d integration steps ", totalTimeSteps);
                
            FPRINTF_CPU;
            FPRINTF_CPU2;
            FPRINTF_CPUD;
            
            fprintf(pFile, "# All quantities except Asy, RG_PE, and sumX are log10(x)\n");
            fprintf(pFile, "# Log of absolute values for E and dE/dt as they can be negative\n");
            fprintf(pFile, "# Units: t and dt in s; E in erg; dE/dt in erg/g/s; others dimensionless \n");
            fprintf(pFile, "#\n");
            
            string str2 = "#      t       dt  1/Rmin   Reaction Rmin    1/Rmax  Reaction Rmax\n";
            fprintf(pFile2, "# All double quantities are log10(x); rates in units of s^-1\n#\n");
            fprintf(pFile2, stringToChar(str2));
            
            for(int i=0; i<plotSteps; i++){
                
                fprintf(pFile2, "%7.4f %7.4f %7.4f %s %7.4f %s\n", 
                    tplot[i], dtplot[i], log10(1.0/slowestRatePlot[i]), 
                    reacLabel[ slowestRateIndexPlot[i]],
                    log10(1.0/fastestRatePlot[i]),
                    reacLabel[ fastestRateIndexPlot[i]]
                );
            }
            
            // Write header for file pointed to by pFile
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[i];
                app.append(Xstring);
                app.append(iso);
                app.append(")     ");
            }
            
            str1.append(app);
            str1.append("\n");
            fprintf(pFile, stringToChar(str1));
            
            fprintf(pFile, "\n");
            
            // Write header for file pointed to by pFile3

            for(int i=0; i<LX; i++){
                iso = isoLabel[i];
                appflux.append(Fpstring);
                appflux.append(iso);
                appflux.append(")   ");
            }
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[i];
                appflux.append(Fmstring);
                appflux.append(iso);
                appflux.append(")   ");
            }
            
            strflux.append(appflux);
            fprintf(pFile3, stringToChar(strflux));
            
            
            // Write the data to the file line by line using concatenated fprintf
            // statements.

            // Loop over timesteps for plot output
            
            for(int i=0; i<plotSteps; i++){
                
                // Initial data fields for t, dt, sumX, fraction of asymptotic
                // isotopes, and fraction of reaction groups in equilibrium.
                
                fprintf(pFile, "%+6.3f %+6.3f %6.3f %6.3f %5.3f %5.3f %5.3f",
                    tplot[i], dtplot[i], EReleasePlot[i], dEReleasePlot[i], 
                    (double)numAsyplot[i]/(double)ISOTOPES,
                    (double)numRG_PEplot[i]/(double)numberRG,
                    sumXplot[i]
                );// to plot diffX add diffXplot[i] and %5.3e at the end of 712
                
                // Now add one data field for each X(i) in plotXlist[]. Add
                // 1e-24 to X in case it is identically zero since we are
                // taking the log.
                
                for(int j=0; j<LX; j++){
                    fprintf(pFile, " %5.3e", log10(Xplot[j][i]+1e-24)); //add log10 in to use log plots with X (set y range in gnufile for gnuplot_X.gnu)
                }
                
                fprintf(pFile, "\n");
                
                // Fluxes
                
                fprintf(pFile3, "\n%+6.3f %+6.3f",// %6.3f %6.3f %5.3f %5.3f %5.3f",
                        tplot[i], dtplot[i]
                );
                
                // Now add one data field for each FplusSumPlot. Add
                // 1e-24 to X in case it is identically zero since we are
                // taking the log.
                
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", log10(abs( FplusSumPlot[j][i]+1e-24) ));
                }
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", log10(abs( FminusSumPlot[j][i]+1e-24)));
                }
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", 
                        log10( abs(FplusSumPlot[j][i] - FminusSumPlot[j][i] + 1e-24) ));
                }
                
            }
            
            // Close output files
            
            fclose (pFile);
            fclose (pFile2);
            fclose (pFile3);
            
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
        
        static int returnNetIndexZN(int z, int n) {
            
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
        // Static function Utilities::isInNet(Z,N) to return true if given (Z,N) 
        // in the network, false otherwise.
        // ----------------------------------------------------------------------
        
        static bool isInNet(int Z, int N) {
            
            if (returnNetIndexZN(Z, N) < 0) {
                return false;
            } else {
                return true;
            }
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::minimumOf(x,y) to return minimum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int minimumOf(int x, int y) {            // integers
            return (x < y) ? x : y;
        }
        
        static double minimumOf(double x, double y) {   // doubles
            return (x < y) ? x : y;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::maximumOf(x,y) to return maximum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int maximumOf(int x, int y) {            // integers
            return (x > y) ? x : y; 
        }
        
        static double maximumOf(double x, double y){    // doubles
            return (x > y) ? x : y; 
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::stringToChar(s) to convert string to char 
        // so that it will print in printf. Assumes max of 50 characters.  
        // Typically a string type can be printed with cout but a string given 
        // to printf usually displays garbage because of type issues in the 
        // C function printf. This function converts a string to a corresponding 
        // character array, which either printf or cout can print. Presently
        // assumes the string has no more than 50 characters.  Change the
        // dimension of cs[] to increase that.
        // ----------------------------------------------------------------------
        
        static char* stringToChar(string s){
            
            static char cs[50];
            strcpy(cs, &s[0]);   // alternatively strcpy(cs, s.c_str());
            return cs;
            
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
        // executing long, pointless loop.
        // ----------------------------------------------------------------------
        
        static void testTimerCPU(){
            
            double a, b;
            
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
    
    // Make data fields private, with external access to them through public setter 
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
        
        // Temperatures in units of 10^9 K for partition function table (see pf[]). 
        const double Tpf[PF] = { 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 
            0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 
            4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f };
    
    public:
        
        // Public setter functions to set values of private class fields
        
        void setisoIndex(int iso){isoindex = iso;}
        
        void setZ(int z){
            ZZ = z;             // Class field
            Z[isoindex] = z;    // Array in main
        }
        
        void setN(int n){
            NN = n;             // Class field
            N[isoindex] = n;    // Array in main
        }
        
        void setA(int a){
            A = a;              // Class field
            AA[isoindex] = a;   // Array in main
        }
        
        void setY(double y){
            YY = y;             // Class field
            XX = (double)A*YY;  // Corresponding X
            Y[isoindex] = y;    // Array in main
            X[isoindex] = XX;   // Array in main
        }
        
        void setY0(double y){ YY0 = y; }
        
        void setX(double x){
            XX = x;             // Class field
            YY = XX/(double)A;  // Corresponding Y
            X[isoindex] = x;    // Array in main
        }
        
        void setM(double m){
            MassExcess = m;             // Class field
            massExcess[isoindex] = m;   // Array in main
        }
        
        void setisoLabel(char *lpt){
            
            // Fill character array.  *lpt points to initial character.
            // All characters of the character array are accessed by
            // incrementing *lpt.
            
            for(int j=0; j<5; j++){
                IsoLabel[j] = *(lpt+j);               // Class field
                isoLabel[isoindex][j] = IsoLabel[j];  // Array in main
            }
            
        }
        
        // Partition function for this species (24 entries)
        
        void setpf(int i, double pfvalue){pf[i] = pfvalue; };
        
        void setfminus(double f){
            fminus = f;                  // Class field
            FminusSum[isoindex] = f;     // Array in main
        }
        
        void setfplus(double f){
            fplus = f;                  // Class field
            FplusSum[isoindex] = f;     // Array in main
        }
        
        void setkeff(double f){ keff=f; }
        
        void setdYdt(double d){
            dYdt = d; 
            dXdt = d*(double)A;
            dYDt[isoindex] = d;
        }
        
        void setdXdt(double d){
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
            printf("\nLabel=%s", object1.getLabel());
          
            char label[5];
            for(k=0; k<5; k++){
               label[k] = object1.getLabel(k);
            }
            printf("\nLabel=%s", label);
            
         assuming that the data fields of object1 have already been populated using
         the setX functions.
        */
        
        char *getLabel() {return IsoLabel; };         // return pointer to label array
        char getLabel(int k) {return IsoLabel[k]; };  // return kth character of array
        
        // Partition function table
        double getpf(int i){return pf[i]; }
        
        // Partition function table temperatures 
        double getTpf(int i){return Tpf[i]; }
        
        double getfminus(){return fminus; }
        
        double getfplus(){return fplus; }
        
        double getkeff(){ return keff; }
        
        double getdYdt(){return dYdt; }
        
        double getdXdt(){return dXdt; }
        
        
};      // End of class Species



/* Class Reaction to describe all reactions in the network.  Create new instance
 * for each reaction in the network. Inherits from Utilities */

class Reaction: public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
    private:
        
        int reacIndex;               // Index of reaction in current run for each (Z,N)
        int reacClass;               // Reaction class for reaclib (1-8)
        int reacGroupClass;          // Reaction group class (1-5 surrogate for A-E)
        int RGmemberIndex;           // index of reaction within its reaction group
        int rgindex;                 // Index of RG (0, 1, ... #RG)
        
        // NOTE: Using the C++ class string instead of char to handle strings 
        // won't compile on my system unless #include <iostream> is included. The
        // C printf command also won't work correctly because it isn't typesafe.
        // See the discussion at
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
        //      printf("\n%s", test+test2);
        //
        // it will likely print garbage. However, you can print a string with printf
        // if it is first converted to a char:
        //
        //      string s = "Howdy World!";
        //      char cs[s.size() + 1];
        //      strcpy(cs, &s[0]);  // or strcpy(cs, s.c_str());
        //      printf("\n\nstring=%s\n", strcpy(cs, &s[0]));
        //
        // The function getreacChar() below returns the string reacString as a
        // pointer to a character array that will work in printf. Alternatively,
        // Utilities::stringToChar() will do same thing.
        
        string reacGroupClassLett;   // Letter equivalent (A-E) for reacGroupClass
        string reacGroupSymbol;      // Schematic equil reaction (e.g. a+b<->c)
        int numberReactants;         // Number species on the left side of reaction
        int numberProducts;          // Number species on the right side of reaction
        int numberIsotopes;          // numberReactants + numberProducts in reaction
        string reacString;           // String describing reaction
        string resonanceType;        // Whether resonant (r) or non-resonant (nr)
        int isEC;                    // Whether electron capture reaction (1) or not (0)
        int isReverse;               // Whether reverse reaction (1) or not (0)
        int ispeforward;             // Whether reactions is labeled "forward" in PE scheme
        bool isEquil;                // Whether in a RG satisfying PE conditions
        
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
        
        // Precomputed temperature factors for ReacLib rates.  Computed in computeTfacs(T9), where
        // T9 is the temperature in units of 10^9 K.
        
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
        double flux;                 // Current net flux of reaction
        double dErate;               // Current rate of energy release
        char cs[20];                 // Utility character array
        char ccs[20];                // Utility character array
        string ss;                   // Utility string
  
  
    public:
        
        // Constructor
        
        Reaction(){
            isEquil = false;
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
                    reacGroupClassLett = "A";
                    reacGroupSymbol = "a<->b";
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 2:
                    reacGroupClassLett = "B";
                    reacGroupSymbol = "a+b<->c";
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 3:
                    reacGroupClassLett = "C";
                    reacGroupSymbol = "a+b+c<->d";
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 4:
                    reacGroupClassLett = "D";
                    reacGroupSymbol = "a+b<->c+d";
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 5:
                    reacGroupClassLett = "E";
                    reacGroupSymbol = "a+b<->c+d+e";
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
            }
        }
        
        void setreacGroupSymbol(string s){reacGroupSymbol = s; }
        
        void setRGmemberIndex(int i){
            
            // Set field of this class
            RGmemberIndex = i;
            
            // Write to corresponding array RMmemberIndex[] in main
            RGMemberIndex[reacIndex] = i;
            
        }
        
        void setrgindex(int i){rgindex = i; }
        
        void setnumberReactants(int i){
            
            numberReactants = i;                 // Field in this class
            NumReactingSpecies[reacIndex] = i;   // Array in main
            
        }
        
        void setnumberProducts(int i){
            
            numberProducts = i;              // Field in this class
            NumProducts[reacIndex] = i;      // Array in main
            
        }
        
        void setreacString(string s){ 
            
            // Set field of this class
            reacString = s; 
            
            // Set corresponding character array reacLabel in main
            char p[s.length()];  
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
                
                // Change array in main
                reacZ[getreacIndex()][i] = z[i];
            }

        }
        
        void setreactantN(int n[]){
            for (int i=0; i<numberReactants; i++){
                
                // Change field in this class (Reaction)
                reactantN[i] = n[i];
                
                // Change array in main
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
                
                // Change array in main
                prodZ[getreacIndex()][i] = z[i];
            }
        }
        
        void setproductN(int n[]){
            for (int i=0; i<numberProducts; i++){
                
                // Change field in this class (Reaction)
                productN[i] = n[i];
                
                // Change array in main
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
                ReactantIndex[reacIndex][i] = n[i];   // Array in main
            }
        }
        
        void setproductIndex(int n[]){
            for (int i=0; i<numberProducts; i++){
                productIndex[i] = n[i];               // Class field
                ProductIndex[reacIndex][i] = n[i];    // Array in main
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes no arguments
        // and constructs isoIndex[] as the concatenation of reactantIndex[]
        // and productIndex[], assuming that those fields have been populated.
        
        void setisoIndex(void){
            for(int i=0; i<numberReactants; i++){
                isoIndex[i] = reactantIndex[i];
            }
            for(int i=0; i<numberProducts; i++){
                isoIndex[i+numberReactants] = productIndex[i];
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes two arguments
        // sets a value of a particular isoIndex[].
        
        void setisoIndex(int i, int j){isoIndex[i] =  j;}
        
        void setdensfac(double d){ densfac = d;}
        
        void setrate(double r){ rate = r; }
        
        void setRrate(double r){ Rrate = r; }
        
        void setflux(double f){ flux = f; }
        
        // Overloaded dErate setter
        
        void setdErate(double r){ dErate = r; }
        void setdErate(void){ dErate = flux*Q; }
        
        
        // Public Reaction getter functions to get values in private fields
        
        int getreacIndex(){ return reacIndex; }

        int getreacClass(){ return reacClass; }
        
        int getreacGroupClass(){ return reacGroupClass; }
        
        // Return reacGroupSymbol as string
        string getreacGroupSymbol(){ return reacGroupSymbol; }
        
        // Return reacGroupSymbol as char array so it will work in printf
        
        char* getreacGroupChar(){
            strcpy(cs, reacGroupSymbol.c_str());
            return cs;
        }
        
        int getRGmemberIndex(){return RGmemberIndex;}
        
        int getrgindex(){return rgindex;}
        
        string getreacGroupClassLett(){ return reacGroupClassLett; }
        
        int getnumberReactants(){ return numberReactants; }
        
        int getnumberProducts(){ return numberProducts; }
        
        // return reacString as string
        string getreacString(){ return reacString; }
        
        // Return reacString as char array so it will work in printf
        char* getreacChar(){
            strcpy(ccs, reacString.c_str());
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
                printf("\n\nERROR: k-1=%d larger than number reactants %d", 
                    k, numberReactants);
                return -1;
            } else {
                return reactantZ[k];
            }
        }
        
        int getreactantN(int k){
            if(k > numberReactants-1){
                printf("\n\nERROR: k-1=%d larger than number reactants %d", 
                    k, numberReactants);
                return -1;
            } else {
                return reactantN[k];
            }
        }
        
        int getreactantA(int i){
           return (reactantZ[i] + reactantN[i]);
        }
        
        int getproductZ(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d", 
                    k, numberProducts);
                return -1;
            } else {
                return productZ[k];
            }
        }
        
        int getproductN(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d", 
                    k, numberProducts);
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
                ss = "\n\nERROR Reaction::getreactantIndex(k): k-1 = %d ";
                ss += "larger than # reactants %d";
                printf(stringToChar(ss), 
                    k, numberReactants);
                return -1;
            } else {
                return reactantIndex[k];
            }
        }
        
        int getproductIndex(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d", 
                    k, numberProducts);
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
            
            // Optional diagnostic output if showFparsing==1
            if(showFparsing == 1){
                fprintf(pFileD, "\n\n--- MAX F+ and F- INDEX FOR EACH ISOTOPE ---\n");  
                for(int i=0; i<numberSpecies; i++){
                    fprintf(pFileD, "\nIsotope index = %d  %s  Max index F+ = %d  Max index F- = %d", 
                        i, isoLabel[i], FplusIsotopeCut[i]-1, FminusIsotopeCut[i]-1);
                }
            }
            
            for(int i=0; i<totalFplus; i++){
                MapFplus[i] = tempInt1[i];
            }
            
            for(int i=0; i<totalFminus; i++){
                MapFminus[i] = tempInt2[i];
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
                for(int j=0; j<SIZE; j++)
                {
                    if(reacMask[i][j] > 0)
                    {
                        FplusFac[tempCountPlus] = (double)reacMask[i][j];
                        tempCountPlus ++;
                    }
                    else if(reacMask[i][j] < 0)
                    {
                        FminusFac[tempCountMinus] = -(double) reacMask[i][j];
                        tempCountMinus ++;
                    }   
                }
            }
            
            // Optional diagnostic output
            
            if(showFparsing == 1){
                fprintf(pFileD, "\n\n\n--- %d NON-VANISHING F+ SOURCE TERMS ---\n", totalFplus);
                fprintf(pFileD, "\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
                       isoLabel[FplusIsotopeIndex[0]], FplusIsotopeIndex[0],
                       numFluxPlus[FplusIsotopeIndex[0]]);
                for(int i=0; i<totalFplus; i++){
                    fprintf(pFileD, "\n   Isotope index = %d F+ index = %d Reac index = %d  %s", 
                           FplusIsotopeIndex[i], i, MapFplus[i], 
                           reacLabel[MapFplus[i]]); 
                    if(i == (FplusIsotopeCut[FplusIsotopeIndex[i]] - 1)  && i != totalFplus-1)
                    {
                        fprintf(pFileD, "\n");
                        fprintf(pFileD, "\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
                               isoLabel[FplusIsotopeIndex[i+1]], FplusIsotopeIndex[i+1],
                               numFluxPlus[FplusIsotopeIndex[i+1]]);
                    }
                }   
                fprintf(pFileD, "\n\n\n--- %d NON-VANISHING F- SOURCE TERMS ---\n", totalFminus);
                fprintf(pFileD, "\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):", 
                       isoLabel[FminusIsotopeIndex[0]], FminusIsotopeIndex[0],
                       numFluxMinus[FminusIsotopeIndex[0]] );
                for(int i=0; i<totalFminus; i++){
                    fprintf(pFileD, "\n   Isotope index = %d F- index = %d Reac index=%d  %s", 
                           FminusIsotopeIndex[i], i, MapFminus[i], reacLabel[MapFminus[i]]);
                    if(i == (FminusIsotopeCut[FminusIsotopeIndex[i]] - 1) && i != totalFminus-1 ){
                        fprintf(pFileD, "\n");
                        fprintf(pFileD, "\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):", 
                               isoLabel[FminusIsotopeIndex[i+1]], FminusIsotopeIndex[i+1],
                               numFluxMinus[FminusIsotopeIndex[i+1]]
                        );
                    }
                }
                fprintf(pFileD, "\n\n");
            }
        }
        
        
        // Function Reaction::populateFplusFminus() to populate F+ and F- for each
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
        // hydro timestep, since the density will generally change over a hydro timestep 
        // in each zone.
        
        void computeConstantFacs(double T9, double rho){

            // Temperature factors in ReacLib rate formula.

            T93 = powf(T9, THIRD); 
            t1 = 1/T9;
            t2 = 1/T93;
            t3 = T93;
            t4 = T9;
            t5 = T93*T93*T93*T93*T93;
            t6 = logf(T9);
            
            // Multiply the statistical prefactor by the appropriate 
            // density factors (1 for 1-body, rho for 2-body, and rho^2 
            // for 3-body reactions).
            
            Dens[0] = 1.0f;
            Dens[1] = rho;
            Dens[2] = rho*rho;
            densfac = prefac * Dens[numberReactants - 1];
            setdensfac(densfac); 
            
        }
        
        // Function to set temperature factors in ReacLib
        // rate formula.  Reset each time T9 changes but
        // stays constant as long as T9 doesn't change.
        
        void computeTfacs(double T9){
            T93 = powf(T9, THIRD); 
            t1 = 1/T9;
            t2 = 1/T93;
            t3 = T93;
            t4 = T9;
            t5 = T93*T93*T93*T93*T93;
            t6 = logf(T9);
        }
        
        // Reaction::computeDensityFactors(double) to multiply the statistical prefactor by 
        // the appropriate density factors (1 for 1-body, rho for 2-body, and rho^2 for 3-body. 
        // This is required at the beginning of each network integration of the hydro timestep, 
        // since the density will generally change over a hydro timestep in each zone.
        
        void computeDensityFactors(double rho){
            
            Dens[0] = 1.0f;
            Dens[1] = rho;
            Dens[2] = rho*rho;
            densfac = prefac * Dens[numberReactants - 1];
            setdensfac(densfac);
        }
        
        // Reaction::computeRate(double, double) to compute rates at T and rho. 
        // The quantity rate is the temperature-dependent part, including a possible 
        // partition-function correction.  The quantity Rrrate is rate multiplied by
        // appropriate density and statistical factors, which give units of s^-1.  The 
        // flux follows from multiplying Rrate by appropriate abundances Y in computeFlux().
        
        void computeRate(double T9, double rho){
            
            // Temperature-dependent rate from ReacLib library
            
            rate = expf( p[0] + t1*p[1] + t2*p[2] + t3*p[3] + t4*p[4] + t5*p[5] + t6*p[6] );
            
            // If necessary, correct rate by multiplying by partition functions
            
            pfUpdate();
            
            setrate(rate);

            // Full rate factor in s^-1 (rate above multiplied by density factors)
            
            Rrate = getdensfac() * rate;
            setRrate(Rrate);
            
            // Write to rate array in main
            
            Rate[getreacIndex()] = Rrate;
            
        }
        
        
        // Function Reaction::pfUpdate() to correct the rates using partition 
        // function factors if appropriate.
        
        void pfUpdate(){
            
            double pfnum;
            double pfden;
            double pfFactor;
            
            // Make a partition function correction if this is reverse reaction in
            // sense defined in ReacLib (defined by field Reaction::isReverse=true). 
            // Realistic calculations at higher temperatures should use
            // the partition functions so generally set dopf=true unless testing.
            // Partition functions are very near 1.000 if T9 < 1, so we will typically
            // only implement partition function correction if T9 > pfCut9 = 1.0, but
            // the table of partition functions allows pfCut9 as small as 0.1.
            // Interpolation is in the log10 of the temperature, so pass log10(T9)
            // rather than T9 to pfInterpolator (index, logt9). Because for the 
            // temperatures of interest the partition functions for all light ions
            // (protons, neutrons, alphas, tritons) are equal to 1.0, the structure
            // of the 8 Reaclib reaction classes specified by Reaction::reacClass
            // means that this correction is only required for reverse reactions
            // in Reaclib classes reacClass = 2, 5.
            
printf("\n********** %s isReverse=%d reacClass=%d", Utilities::stringToChar(reacString), isReverse, reacClass );
            
            if(dopf && T9 > pfCut9 && isReverse){
                
                if(reacClass == 2){
                    
                    pfden = pfInterpolator (reactantIndex[0], log10(T9));
                    pfnum = pfInterpolator (productIndex[1], log10(T9));
                    
                } else if(reacClass == 5){
                    
                    pfden = pfInterpolator (reactantIndex[1], log10(T9));
                    pfnum = pfInterpolator (productIndex[1], log10(T9));
                    
                } else {
                    
                    pfden = 1.0;
                    pfnum = 1.0;
                    
                }
                
                pfFactor = pfnum/pfden;
                rate *= pfFactor;
                
printf("\n           pfnum=%7.4e pfden=%7.4e pfFactor=%7.4e newrate=%7.4e oldrate=%7.4e", pfnum, pfden, pfFactor, rate, rate/pfFactor);
                
            }
        }
        
        
        // ------------------------------------------------------------------------
        // Return partition function of isotope labeled by isoIndex at log_10 of
        // temperature T9. Note that the 2nd argument is log10(T9), not T9,
        // because the interpolation in the partition function table is in the 
        // log10 of the temperature.  The following commented-out code assumes
        // that the object interpolatepf of the SplineInterpolator class has
        // first invoked the interpolatepf.bisection method to use bisection 
        // to find the interval containing root and store the lower index of
        // that interval in lowPFindex. Then SplineInterpolator interpolates
        // the root restricted to that interval.  This guards against the
        // spline interpolator finding the wrong root if there are multiple
        // roots (as could be true in the general case, though probably not here
        // since the function is typically monotonic).
        // ------------------------------------------------------------------------
        
        double pfInterpolator(int index, double logt9) {
            
            // Following commented out for testing purposes until spline interpolator
            // is implemented
            
//             double rdt;
//             double term1;
//             double term2;
//             double sumterms;
//             double bob;
//             rdt = (logt9 - Tpf[lowPFindex]) / (Tpf[lowPFindex + 1] - Tpf[lowPFindex]);
//             term1 = rdt * Math.log(pf[Z][N][lowPFindex + 1]);
//             term2 = (1.0 - rdt) * Math.log(pf[Z][N][lowPFindex]);
//             sumterms = term1 + term2;
//             bob = Math.exp(sumterms);
//             // System.out.println("PF stuff: "+t9+" "+Z+" "+N+" "+rdt+" "+sumterms+" "+bob);
//             return bob;
            
            return 1.0;  // Temporary
            
        }
        
        // Function Reaction::showRates() to display computed rates for this
        // Reaction object.
        
        void showRates(){
            
            fprintf(pFileD, "\n%d %19s RG=%d densfac=%6.3e rate= %8.5e Rrate=%8.5e", 
                   getreacIndex(), getreacChar(), getreacGroupClass(), getdensfac(), 
                   getrate(), getRrate()
            );
            
        }
        
        
        // Function Reaction::computeFlux() to compute the current flux for reaction 
        // corresponding to this Reaction object.
        
        void computeFlux(){
            
            // If we are imposing partial equilibrium and this reaction is part of a 
            // reaction group that is in partial equilibrium, set its flux to zero 
            // and return.
            
            if(doPE  && t > equilibrateTime && !reacIsActive[reacIndex]){
                flux = 0.0;
                Flux[reacIndex] = flux;     // Put in flux array in main
                return;
            }
            
            // Otherwise, this reaction is in a RG that is not in equilibrium, so we
            // need to compute its flux
            
            string s;
            double kfac;
            
            switch(numberReactants){
                
                case 1:    // 1-body reactions
                    
                    kfac = Rrate;
                    flux = kfac*Y[ reactantIndex[0] ];  
                    Flux[reacIndex] = flux;         // Put in flux array in main
                    fastSlowRates(kfac);
                    
                    if(showFluxCalc == 1){
                        fprintf(pFileD, "\n%d %18s reactants=%d iso0=%d Rrate=%7.3e Y1=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0],  Rrate, 
                            Y[ reactantIndex[0] ], flux);
                    }
                    
                    break;
                    
                case 2:    // 2-body reactions  
                    
                    kfac = Rrate * Y[ reactantIndex[0] ];
                    flux = kfac * Y[ reactantIndex[1] ];    
                    Flux[reacIndex] = flux;         // Put in flux array in main
                    fastSlowRates(kfac);
                    
                    if(showFluxCalc == 1){
                        fprintf(pFileD, "\n%d %18s reactants=%d iso0=%d iso1=%d Rrate=%7.3e Y1=%7.3e Y2=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], 
                            reactantIndex[1], Rrate, Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], flux);
                    }
                    
                    break;
                    
                case 3:    // 3-body reactions
                    
                    kfac = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ];
                    flux = kfac * Y[ reactantIndex[2] ];
                    Flux[reacIndex] = flux;         // Put in flux array in main
                    fastSlowRates(kfac);
                    
                    if(showFluxCalc == 1){
                        fprintf(pFileD, 
                            "\n%d %18s reactants=%d iso0=%d iso1=%d iso2=%d Rrate=%7.3e Y1=%7.3e Y2=%7.3e Y3=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], reactantIndex[1], 
                            reactantIndex[2], Rrate, Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], 
                            Y[ reactantIndex[2] ], flux);
                    }
                    
                    break;
            }
            
        }    // End of function computeFlux()
        
      
      
        // Reaction::fastSlowRates(double) to store fastest and slowest rates 
        // in this timestep, and overall in the calculation. These rates are the
        // rates kfac computed in computeFlux().
        
        void fastSlowRates(double testRate){
            
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
            
        }  // End of function fastSlowRates()
        
        
        // Function Reaction::sumFplusFminus() to sum the total F+ and F- for each isotope.  
        
        static void sumFplusFminus(){
            
            int minny = 0;
            double accum;
            double dydt;
            
            for(int i=0; i < numberSpecies; i++){   
                
                // Sum F+ for each isotope
                
                if(i>0) minny = FplusMax[i-1]+1;
                if(showFluxCalc == 1) fprintf(pFileD, "\n\nFplusMax=%d", FplusMax[i-1]);
                accum = 0.0;    
                for(int j=minny; j<=FplusMax[i]; j++){
                    accum += Fplus[j];
                    if(showFluxCalc == 1) fprintf(pFileD, "\ni=%d j=%d Fplus=%g FplusSum=%8.4e", 
                        i, j, Fplus[j], accum);
                }
                
                setSpeciesfplus(i, accum);        // Also sets FplusSum[i] = accum;
                
                // Sum F- for each isotope
                
                minny = 0;
                if(i>0) minny = FminusMax[i-1]+1;
                if(showFluxCalc == 1) fprintf(pFileD, "\n\nFminusMax=%d", FminusMax[i-1]);
                accum = 0.0;
                for(int j=minny; j<=FminusMax[i]; j++){
                    accum += Fminus[j];
                    if(showFluxCalc == 1) fprintf(pFileD, "\ni=%d j=%d Fminus=%g FminusSum=%8.4e", 
                        i, j, Fminus[j], accum);
                }
                
                setSpeciesfminus(i, accum);      // Also sets FminusSum[i] = accum and keff
                setSpeciesdYdt(i, FplusSum[i] - FminusSum[i]);
                
            }
            
        }
                
};  // End class Reaction



// Class ReactionVector with  functions that create reaction vectors. Inherits from 
// class Utilities.


class ReactionVector:  public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
    private:
        
    public:
        
        /* Use the information constructed in parseF (in particular reacMask[j][k]) to create the 
        reaction vectors for the network using the static function 
        ReactionVector::makeReactionVectors(). */
        
        static void makeReactionVectors(){
            
            // Write out the array containing components of the reaction vectors
            
            int uppity = minimumOf(25, numberSpecies);  // limit printout width to 25 species
            fprintf(pFileD, "\n\nREACTION VECTOR ARRAY (%d Reaction vectors with %d species components):\n", 
                numberReactions, ISOTOPES);
            fprintf(pFileD, "\nReaction \\ Species           ");
            for(int k=0; k<uppity; k++){
                fprintf(pFileD, "%4s  ",isoLabel[k]);
            }
            fprintf(pFileD, "\n");
            for(int j=0; j<numberReactions; j++){
                fprintf(pFileD, "%4d %22s [ ",j,reacLabel[j]);
                for(int k=0; k<uppity-1; k++){
                    fprintf(pFileD, "%2d    ", reacMask[k][j]);
                }
                fprintf(pFileD, "%2d ]\n", reacMask[uppity-1][j]);
            }
            
            // -----------------------------------------------------------------------
            // Now implement reaction vectors as GSL vectors so that we can use the
            // GSL and GSL_BLAS API to manipulate them.  The prototypes for doing this
            // may be found in the example code arrayPointers.c and matrixGSL.c
            // -----------------------------------------------------------------------
            
            rvPt = rv;       // Set pointer to beginning address of array rv
            
            // Allocate an array populated with GSL vectors
            
            if(showRVdetails == 1) printf("\nAllocating an array rv[] of GSL vectors\n");
            
            for(int i=0; i<SIZE; i++){
                
                //Prototype GSL reaction vector
                
                gsl_vector * v1 = gsl_vector_alloc (ISOTOPES); 
                
                // Set elements of rv[] pointed to by *rvPt equal to GSL vectors
                
                *(rvPt+i) = *v1;   
            }
            
            // Fill vector component entries created above with data contained in  
            // reacMask[j][i] (notice reversed indices)
            
            if(showRVdetails == 1){
                fprintf(pFileD, "\nPopulate vector components of array rv[i]_j with reacMask[j][i]:");
            }
            
            for (int i = 0; i < SIZE; i++) {
                if(showRVdetails == 1) printf("\n\nrv[%d]",i);
                for(int j=0; j<ISOTOPES; j++){
                    gsl_vector_set (rvPt+i, j, reacMask[j][i]);
                    
                    // Retrieve the vector component just stored and print it
                    
                    int temp = gsl_vector_get(rvPt+i, j);
                    if(showRVdetails == 1){
                        fprintf(pFileD, "\ni=%d j=%d  reacMask[%d][%d] =%3d  rv[%d]_%d =%3d", 
                            i, j, i, j, reacMask[j][i], i, j, temp);
                    }
                }
            }
        
            
            // Display reaction vectors as component list
            
            fprintf(pFileD, "\nGSL REACTION VECTOR COMPONENTS (%d reaction vectors with %d components)\n",
                SIZE, ISOTOPES);
            
            for (int i = 0; i < SIZE; i++) {
                fprintf(pFileD, "\nrv[%d]: [",i);
                for(int j=0; j<ISOTOPES; j++){
                    
                    // Define a pointer that will point to the GSL vector in array entry rv[i]
                    
                    gsl_vector *vector = rvPt+i;
                    
                    // Assign the jth component of the vector in rv[i] to a variable
                    
                    int component = gsl_vector_get (vector, j);
                    fprintf (pFileD, "%3d", component);
                }
                
                fprintf(pFileD, " ]");
            }
            
        }  // End function makeReactionVectors()
    
    
        
        // ------------------------------------------------------------------------
        // ReactionVector::compareGSLvectors(rv1, rv2) to compare two GSL vectors 
        // of same length.  Returns 0 if they are not equivalent, 1 if they are 
        // the same, 2 if one vector is the negative of the other. The arguments 
        // are pointers to the two GSL vectors.
        // ------------------------------------------------------------------------
        
        int static compareGSLvectors(gsl_vector* rv1, gsl_vector* rv2){
            
            int k, kk;
            
            // Compare rv1 and rv2. Function gsl_vector_equal(rv1, rv2) returns 1 
            // if vectors are equal and 0 if they are not.
            
            k = gsl_vector_equal(rv1, rv2);
            
            if (k==1) return k;  // rv1 = rv2, so return with 1
            
            // Since rv1 and rv2 are not equal, compare rv1  and -rv2
            
            gsl_vector * rv2minus = gsl_vector_alloc(ISOTOPES);
            gsl_vector_memcpy(rv2minus, rv2);
            gsl_vector_scale(rv2minus, -1);
            kk = gsl_vector_equal(rv1, rv2minus);
            
            if(kk==0){
                return 0;  // rv1 not equal to rv2 and not equal to -rv2
            } else {
                return 2;  // rv1 equal to -rv2
            }
            
        }    // End function compareGSLvectors
    
    
        
        // ------------------------------------------------------------------------
        // ReactionVector::sortReactionGroups() uses compareGSLvectors to sort all 
        // reactions in the network into reaction groups labeled by a series of 
        // integers 0, 1, ...  All reactions in a reaction group have the same 
        // reaction vector up to a sign. The array RGindex[] of dimension SIZE 
        // holds the integer labeling reaction group (0, 1, ... #RG) for each 
        // reaction after this function is executed. 
        // ------------------------------------------------------------------------
        
        static void sortReactionGroups(void){
            
            // Cycle over all reaction vectors and compare them pairwise to 
            // assign to reaction groups. The pointer rvPt points to the array
            // rv[] of GSL reaction vectors.
            
            // The integer rindex labels the reaction group
            
            int rindex = -1;
            int ck = -1;
            
            // Initialize
            
            for(int i=0; i<SIZE; i++){
                RGindex[i] = -1;
            }
            
            if(showRGsorting == 1) fprintf(pFileD, "\n\n\n--- SORTING REACTION GROUPS ---");
            
            int scorekeeper;
            for (int i=0; i<SIZE; i++){
                scorekeeper = 0;
                if(i==0) rindex ++;
                if(showRGsorting == 1) fprintf(pFileD, "\n\nRG=%d", rindex);
                for(int j=0; j<SIZE; j++){
                    
                    if(RGindex[i] < 0) RGindex[i] = rindex;
                    ck = compareGSLvectors(rvPt+i, rvPt+j);
                    
                    if(ck > 0 && RGindex[j]< 0) {
                        RGindex[j] = rindex;
                        scorekeeper ++;
                    }
                    if(showRGsorting==1){
                        fprintf(pFileD, "\ni=%d %s j=%d %s RGindex[%d]=%d ck=%d rindex=%d scorekeeper=%d", 
                            i, reacLabel[i], j, reacLabel[j], j, RGindex[j], ck, rindex, scorekeeper);
                    }
                }
                
                // If scorekeeper > 0, this is a reaction group with scorekeeper+1 members, all having
                // the same reaction vector, up to a sign.
                
                if(scorekeeper > 0){
                    
                    // Store the number of reactions in this reaction group for later use
                    
                    RGnumberMembers[rindex] = scorekeeper+1;
                    
                    // Increment the RG number
                    
                    rindex++;
                    
                    if(showRGsorting==1) fprintf(pFileD, "\nFound RG=%d RGnumberMembers=%d", 
                        rindex-1, RGnumberMembers[rindex-1]);
                }

            }
            
            numberRG = rindex;   // Store total number of reaction groups
            
            // Diagnostic showing reaction group associated with each reaction
            
            if(showRGsorting == 1){
                fprintf(pFileD, "\n\n-- SUMMARY OF REACTION GROUPS:\n");
                for(int i=0; i<SIZE; i++){
                    fprintf(pFileD, "\nreaction=%d  %18s RGindex=%d RGmemberIndex=%d", 
                        i, reacLabel[i], RGindex[i], RGMemberIndex[i]);
                }
            }
            
            // Write out the components of the reaction groups
            
            fprintf(pFileD, "\n\n\nPARTIAL EQUILIBRIUM REACTION GROUPS");
            for(int i=0; i<numberRG; i++){
                fprintf(pFileD, "\n\nReaction Group %d:", i);
                int rgindex = -1;
                for(int j=0; j<SIZE; j++){
                    if(RGindex[j] == i){
                        rgindex ++; 
                        setRG(j, RGclass[j], RGindex[j]);
                        fprintf(pFileD, "\n%s reacIndex=%d RGindex=%d RG=%d RGreacIndex=%d isForward=%d RG: %s", 
                            reacLabel[j], j, rgindex, RGclass[j], RGMemberIndex[j],
                            isPEforward[j], stringToChar(RGstring[j]));
                    }
                }
            }
            
            fprintf(pFileD, "\n");
            
        }      // End function sortReactionGroups()
        
    
    
        /*
        * ReactionVector::parseF() to find contributions to F+ and F- of each reaction for each 
        * isotope. This is executed only once at the beginning of the entire calculation to determine 
        * the structure of the network. Since executed only once, make it static so it can
        * be called directly from the class using ReactionVector::parseF(), without having
        * to instantiate.
        */
        
        static void parseF(){
            
            if(showParsing == 1)
                fprintf(pFileD, "\n\n--- Use parseF() to find F+ and F- flux components for each species ---");
            
            int incrementPlus = 0;
            int incrementMinus = 0;
            
            // Loop over all isotopes in the network
            
            for(int i=0; i<numberSpecies; i++) {
                int total = 0;
                int numFplus = 0;
                int numFminus = 0;
                if(showParsing == 1) fprintf(pFileD, "\n");
                
                // Loop over all possible reactions for this isotope, finding those that
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
                        reacMask[i][j] = -total;
                        tempInt2[incrementMinus + numFminus-1] = j;
                        if(showParsing == 1)
                            fprintf(pFileD, "\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F-", 
                                isoLabel[i], j, reacLabel[j], NumReactingSpecies[j], NumProducts[j], totalL, 
                                totalR, total);
                    } 
                    else if(total < 0){          // Contributes to F+ for this isotope
                        numFplus ++;
                        reacMask[i][j] = -total;
                        tempInt1[incrementPlus + numFplus-1] = j;
                        if(showParsing == 1)
                            fprintf(pFileD, "\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F+", 
                                isoLabel[i], j, reacLabel[j], NumReactingSpecies[j], NumProducts[j], totalL, 
                                totalR, total);
                    } else {           // Does not contribute to flux for this isotope
                        reacMask[i][j] = 0;
                    }
                }
                
                // Keep track of the total number of F+ and F- terms in the network for all isotopes
                
                totalFplus += numFplus;
                totalFminus += numFminus;
                
                numFluxPlus[i] = numFplus;
                numFluxMinus[i] = numFminus;
                
                incrementPlus += numFplus;
                incrementMinus += numFminus;
                
                if(showParsing == 1)
                fprintf(pFileD, "\nSpecies=%d %s numF+ = %d numF- = %d", i, isoLabel[i], numFplus, numFminus);
            }
            
            // Display isotope component array
            
            fprintf(pFileD, 
                "\n\n\nFLUX-ISOTOPE COMPONENT ARRAY (negative n for F-; positive n for F+ for given isotope):");
            fprintf(pFileD, "\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
            
            int uppity = minimumOf(30, numberSpecies);  // limit printout width to 30 species
            fprintf(pFileD, "\n\nIndex             Reaction");
            for(int k=0; k<uppity; k++){
                fprintf(pFileD, "%5s", isoLabel[k]);
            }
            for(int j=0; j<numberReactions; j++){
                fprintf(pFileD, "\n%3d %22s",j,reacLabel[j]);
                for(int k=0; k<uppity; k++){
                    fprintf(pFileD, " %4d",reacMask[k][j]);
                }
            }
            
            fprintf(pFileD, 
                "\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.\n", 
                totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
            
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
                gsl_vector_set(abundances, i, a[i]);
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

        // Free allocated matrix and vector
        
        void freeGSL(){
            gsl_vector_free(abundances);
            gsl_matrix_free(fluxes);
        }

};  // end of class MatrixUtils



// Class ReactionGroup to handle reaction groups.  Inherits from class Utilities

class ReactionGroup:  public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions.  Its static functions can be called directly from the class
    // without having to instantiate.
    
    private:
        
        static const int maxreac = 10;         // Max possible reactions in this RG instance
        int nspecies[5] = { 2, 3, 4, 4, 5 };   // Number isotopic species in 5 RG classes
        int niso;                              // Number isotopic species in this RG object
        int RGn;                               // Index this object in RG array (0,1,... #RG)
        int numberMemberReactions;             // Number of reactions in this RG instance
        int memberReactions[maxreac];          // reacIndex of reactions in reaction group
        int numberReactants[maxreac];          // Number of reactants for each reaction in RG
        int numberProducts[maxreac];           // Number of products for each reaction in RG
        int refreac = -1;                      // Ref. reaction for this RG in memberReactions

        int rgclass;                           // Reaction group class (1-5)
        bool isEquil;                          // True if RG in equilibrium; false otherwise
        bool isForward[maxreac];               // Whether reaction in RG labeled forward
        double flux[maxreac];                  // Current flux for each reaction in RG
        double netflux;                        // Net flux for the entire reaction group
        char reaclabel[maxreac][LABELSIZE];    // Member reaction label
        
        // Partial equilibrium quantities
        
        double crg[4];                 // Constants c1,c2, ... (allocate dynamically?)
        int numberC;                   // Number constants crg[] for this RG object
        double rgkf;                   // Forward rate parameter for partial equilibrium
        double rgkr;                   // Reverse rate parameter for partial equilibrium
        
        double aa, bb, cc;             // Quadratic coefficients a, b, c
        double alpha, beta, gamma;     // Coefficients for cubic ~ quadratic approximation
        double qq;                     // q = 4ac-b^2
        double rootq;                  // Math.sqrt(-q)
        double tau;                    // Timescale for equilibrium

        double equilRatio;             // Equilibrium ratio of abundances
        double kratio;                 // Ratio k_r/k_f. Equal to equilRatio at equilibrium
        double eqcheck[5];             // Population ratio to check equilibrium

        double lambda;                 // Progress variable for reaction pair
        double lambdaEq;               // Equilibrium value of progress variable
        
        int reactantIsoIndex[3];       // Species index of reactants
        int productIsoIndex[4];        // Species index of products
        
        int isoindex[5];               // Species index for participants in reaction   
        char isolabel[5][5];           // Isotopic label of species in RG reactions
        int isoZ[5];                   // Z for niso isotopes in the reactions of the group
        int isoN[5];                   // N for niso isotopes in the reactions of the group
        double isoA[5];                // A for niso isotopes in the reactions of the group
        double isoYeq[5];              // Y_eq for niso isotopes in reactions of the group
        double isoY[5];                // Current Y for niso isotopes in reactions of group
        double isoY0[5];               // Y0 for niso isotopes in the reactions of the group

    
    public:
    
    // Constructor
        
    ReactionGroup(int rgn){
        RGn = rgn;
        isEquil = false;
    }
    
    // Public ReactionGroup setter functions to set values of private class fields
    
    void setnumberMemberReactions(int n){numberMemberReactions = n;}
    
    void setmemberReactions (int i, int index){memberReactions[i] = index;}
    
    void setnumberReactants(int i, int j){numberReactants[i] = j;}
    
    void setnumberProducts(int i, int j){numberProducts[i] = j;}
    
    void setisEquil(bool b) {isEquil = b;}
    
    void setisForward(int i, bool b){isForward[i] = b;}
    
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
        // all the forward reactions in a reaction group, for example).
        
        if (refreac == -1) {
            refreac = 0;
            fprintf(pFileD, "\n*** Reaction group %d has no forward reactions ***", 
                RGn);
        }
        
        return refreac;
    }
    
    void setflux(int i, double f){flux[i] = f;}
    
    void setRGn(int rgn){RGn = rgn;}
    
    void setrgclass(int rc){rgclass = rc;}

    void setreaclabel(int k, string s){ 
        
        // Convert from string to char array
        
        char p[s.length()+1];  
        for (int i = 0; i < sizeof(p); i++) { 
            p[i] = s[i]; 
            reaclabel[k][i] = p[i];
        }
    }
       
    // ReactionGroup function to set all fluxes in RG from the array Flux[].
    
    void setRGfluxes(){

        for(int i=0; i<numberMemberReactions; i++){
            setflux(i, Flux[ memberReactions[i] ]);
        }
    }
    
    void setnetflux(double f){netflux = f;}
    
    void setniso(int rgclass){niso = nspecies[rgclass-1];}
    
    void setisoindex(int i, int j){isoindex[i] = j;}
    
    void setisolabel(int k, char s[]){ 
        for (int i = 0; i < 5; i++) { 
            isolabel[k][i] = s[i];
        }
    }
    
    void setisoZ(int i, int j){isoZ[i] = j;}
    
    void setisoN(int i, int j){isoN[i] = j;}
    
    void setisoA(int i, int j){isoA[i] = j;}
    
    void setisoY0(int i, double y){isoY0[i] = y;}
    
    void setisoY(int i, double y){isoY[i] = y;}
    
    void setisoYeq(int i, double y){isoYeq[i] = y;}
    
    void setreactantIsoIndex(int i, int j){
        reactantIsoIndex[i] = j;
    }

    void setproductIsoIndex(int i, int j){
        productIsoIndex[i] = j;
    }
    
    void seteqcheck(int k, double eq){eqcheck[k] = eq;}
    
    
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
    
    char* getreacString(int k){return reaclabel[k];}
    
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
    
    
    // Function ReactionGroup::showRGfluxes to show all current fluxes in 
    // reaction groups
    
    void showRGfluxes(){
        
        fprintf(pFileD, "\n\nRG=%d", RGn);
        
        double fac;
        for(int i=0; i<numberMemberReactions; i++){
            if(isForward[i]){
                fac = 1.0;
            } else {
                fac = -1.0;
            }
            fprintf(pFileD, "\nshowRGfluxes: %d %s RGclass=%d isForward=%d t=%7.4e dt=%7.4e flux=%7.4e", 
                i, reacLabel[memberReactions[i]], rgclass, isForward[i], t, dt, 
                fac*flux[i]
            );
        }
        
        fprintf(pFileD, "\n");
        if(isEquil){
            fprintf(pFileD, "showRGfluxes: NetRGflux=%7.4e\nEQUILIBRATED",  netflux); 
        } else {
            fprintf(pFileD, "showRGfluxes: NetRGflux=%7.4e\nNOT EQUILIBRATED", netflux); 
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
        
        // Compute net flux and progress variable
        
        if (isEquil) {
            netflux = sumRGfluxes();
            lambda = netflux * deltaTime;
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
    // Function ReactionGroup::putY0() to put the values of Y0 at beginning of 
    // timestep into the Y0[] array for this object
    // -----------------------------------------------------------------------
    
    void putY0() {
        
        int ii;
        if(diagnose2) fprintf(pFileD, "\n");
        for (int k = 0; k < niso; k++) {
            ii = isoindex[k];
            isoY0[k] = Y0[ii];
            isoY[k] = isoY0[k];
            if(diagnose2)
            fprintf(pFileD, "\nputY0: t=%8.5e RG=%d niso=%d k=%d isoindex=%d isoY0[%s]=%8.5e isoY0=%8.5e Y[ii]=%8.5e", 
                t, RGn, niso, k, ii, isoLabel[ii],  isoY[k], isoY0[k], Y[ii]);
        }
        if(diagnose2) fprintf(pFileD, "\n");
    }
    
   
   // -----------------------------------------------------------------------
   // Function ReactionGroup::computeC() to compute values of constants crg[]
   // -----------------------------------------------------------------------
   
    void computeC() {
        
        if(diagnose2) fprintf(pFileD, "\n\nRG=%d", RGn);
        
        switch (rgclass) {
            
            // Reaclib class 7, which can't equilibrate because
            // no inverse reactions in library
            
            case -1:   
                break;
                
            case 1:    // a <-> b
                
                crg[0] = isoY0[0] + isoY0[1];
                break;
                
            case 2:    // a+b <-> c
                
                crg[0] = isoY0[1] - isoY0[0];
                crg[1] = isoY0[1] + isoY0[2];
                
                if(diagnose2){
                    
                    fprintf(pFileD, "\ncomputeC: t=%7.4e RG=%d isoY0[0]=%7.4e isoY0[1]=%7.4e isoY0[2]=%7.4e",
                           t, RGn, isoY0[0], isoY0[1], isoY0[2]
                    );
                    fprintf(pFileD, "\ncomputeC: t=%7.4e RG=%d crg[0]=%7.4e crg[1]=%7.4e",
                           t, RGn, crg[0], crg[1]
                    );
                }
                
                break;
                
            case 3:    // a+b+c <-> d
                
                crg[0] = isoY0[0] - isoY0[1];
                crg[1] = isoY0[0] - isoY0[2];
                crg[2] = THIRD * (isoY0[0] + isoY0[1] + isoY0[2]) + isoY0[3];
                
                if(diagnose2){
                    fprintf(pFileD, "\ncomputeC: t=%7.4e RG=%d isoY0[0]=%8.5e isoY0[1]=%8.5e isoY0[2]=%8.5e isoY0[3]=%8.5e",
                            t, RGn, isoY0[0], isoY0[1], isoY0[2], isoY[3]
                    );
                    fprintf(pFileD, "\ncomputeC: t=%7.4e RG=%d crg[0]=%8.5e crg[1]=%8.5e crg[2]=%8.5e",
                            t, RGn, crg[0], crg[1], crg[2]
                    );
                }
                
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
            
            // Reaclib class 7, which can't equilibrate in ReacLib because 
            // no inverse in reaction library
            
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
                
                if(diagnose2)
                fprintf(pFileD, "\ncomputeQuad: t=%7.4e RG=%d aa=%7.4e bb=%7.4e cc=%7.4e",
                    t, RGn, aa, bb, cc);
                       
                break;
                
            case 3:  // a+b+c <-> d
                
                aa = -rgkf * isoY0[0] + rgkf * (crg[0] + crg[1]);
                bb = -(rgkf * crg[0] * crg[1] + rgkr);
                cc = rgkr * (crg[2] + THIRD * (crg[0] + crg[1]));
                
                if(diagnose2)
                fprintf(pFileD, "\ncomputeQuad: t=%7.4e RG=%d aa=%7.4e bb=%7.4e cc=%7.4e",
                    t, RGn, aa, bb, cc);
                
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
        
        // Compute the q = 4ac - b^2 parameter, equil timescale tau, and
        // isoYeq[0] (which is then be used to compute the other isoYeq[].
        
        if (rgclass > 1) {
            qq = computeq(aa, bb, cc);
            rootq = sqrt(-qq + 1.0e-30);
            if (numberMemberReactions > 1) {
                tau = 1 / rootq;
            }
            isoYeq[0] = computeYeq(aa, bb, rootq);
            
            if(diagnose2)
            fprintf(pFileD, "\ncomputeQuad: t=%7.4e RG=%d qq=%7.4e tau=%7.4e isoYeq[0]=%7.4e",
                t, RGn, qq, tau, isoYeq[0]);
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
                
                if(diagnose2)
                fprintf(pFileD, "\ncomputeQuad: t=%6.4e RG=%d isoYeq[0]=%6.4e isoYeq[1]=%6.4e isoYeq[2]=%6.4e eqRatio=%5.3e",
                    t, RGn, isoYeq[0], isoYeq[1], isoYeq[2], equilRatio);

                break;
                
            case 3:    // a+b+c <-> d
                
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = isoYeq[0] - crg[1];
                isoYeq[3] = crg[2] - isoYeq[0] + THIRD * (crg[0] + crg[1]);
                equilRatio = isoY[0] * isoY[1] * isoY[2] / isoY[3];
                
                if(diagnose2)
                fprintf(pFileD, "\ncomputeQuad: t=%5.3e RG=%d isoYeq[0]=%5.3e isoYeq[1]=%5.3e isoYeq[2]=%5.3e isoYeq[3]=%5.3e eqRatio=%5.3e",
                    t, RGn, isoYeq[0], isoYeq[1], isoYeq[2], isoYeq[3], equilRatio);
                
                break;
                
            case 4:    // a+b <-> c+d
                
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = crg[1] - isoYeq[0];
                isoYeq[3] = crg[2] - isoYeq[0];
                equilRatio = isoY[0] * isoY[1] / (isoY[2] * isoY[3]);
                break;
                
            case 5:    //  a+b <-> c+d+e
                
                isoYeq[1] = isoYeq[0] - crg[1];
                isoYeq[2] = alpha - isoYeq[0];
                isoYeq[3] = beta - isoYeq[0];
                isoYeq[4] = gamma - isoYeq[0];
                equilRatio = isoY[0] * isoY[1] / (isoY[2] * isoY[3] * isoY[4]);
                break;
        }
        
        
        // Compute the equilibrium value of the progress variable
        
        lambdaEq = isoY0[0] - isoYeq[0];
        
        // Compute the population ratios used to check equilibration
        
        computeEqRatios();
        kratio = rgkr / rgkf;
        
    }    // End function computeQuad()
    
    
    // ---------------------------------------------------------------------
    // Function ReactionGroup::computeq to compute q = 4ac-b^2 for quadratic 
    // solution
    // ---------------------------------------------------------------------
    
    double computeq(double a, double b, double c) {
        return 4 * a * c - b * b;
    }
    
    // ---------------------------------------------------------------------
    //  Function ReactionGroup::computeYeq to compute Yeq[0]
    // ---------------------------------------------------------------------
    
    double computeYeq(double a, double b, double rootq) {
        return -0.5 * (b + rootq) / a;
    }
    
    
    // ---------------------------------------------------------------------
    // Function ReactionGroup::computeEqRatios to compute array of population 
    // ratios used to check equilibration
    // ---------------------------------------------------------------------
    
    void computeEqRatios() {
        
        // Add 1e-24 to denominator in following to prevent possible divide by zero
        
        double thisDevious = abs((equilRatio - kratio) / (kratio + 1.0e-24));
        
        if (isEquil && thisDevious > mostDevious) {
            mostDevious = thisDevious;
            mostDeviousIndex = RGn;
        }
        
        if(diagnose2)
        fprintf(pFileD, "\ncomputeEqRatios: t=%6.4e RG=%d equilRatio=%6.4e kratio=%6.4e thisDev=%6.4e mostDev=%7.4e equil=%d",
            t, RGn, equilRatio, kratio, thisDevious, mostDevious, isEquil);
        
        // The return statements in the following if-clauses cause reaction
        // groups already in equilibrium to stay in equilibrium. If the 
        // deviousMax > tolerance check is implemented it can cause a
        // reaction group to drop out of equilibrium.
        
        if (isEquil && thisDevious < deviousMax) {
            return;
        } else if (isEquil && thisDevious > deviousMax && doPE
            && t > equilibrateTime) {
            removeFromEquilibrium();
            return;
        }
            
        Yminner = 1000;
        maxeqcheck = 0;
        mineqcheck = 1000;
        
        // Determine if reaction group RG is in equilibrium: set isEquil to default value
        // of true and then try to falsify
        
        isEquil = true;
        
        for (int i = 0; i < niso; i++) {
            
            // Note: something like the following probably required because
            // otherwise we will divide by zero for isotopes early in the 
            // calculation that have no population.
            
            if (isoYeq[i] == 0 || isoY[i] == 0) {
                isEquil = false;
                break;
            }
            
            eqcheck[i] = abs(isoY[i] - isoYeq[i]) / isoYeq[i];
            
            if(t > equilibrateTime && diagnose2) 
            fprintf(pFileD, 
                "\ncomputeEqRatios: iso=%d %s RG=%d t=%7.4e isoYeq=%7.4e isoY=%7.4e eqcheck=%7.4e",
                i, isolabel[i], RGn, t, isoYeq[i], isoY[i], eqcheck[i] );
            
            // Store some min and max values
            
            if (eqcheck[i] < mineqcheck)
                mineqcheck = eqcheck[i];
            if (eqcheck[i] > maxeqcheck)
                maxeqcheck = eqcheck[i];
            if (isoYeq[i] < Yminner)
                Yminner = isoYeq[i];
            
            // Set equilibrium to false if any eqcheck[] greater than
            // tolerance, or if equilibrium abundance is small relative 
            // to equiTol (which can cause numerical issues in judging 
            // whether in equilibrium).
            
            if (t < equilibrateTime || eqcheck[i] > equiTol || isoYeq[i] < Ythresh) {
                
                isEquil = false;
            
                // break; // Note: this break won't affect results, but
                // would affect diagnostic values of eqcheck[]
                // since they will all be zero after the break.
            
            }
        }
        
        // Set isEquil field of network species vectors to true if isotope
        // participating in equilibrium
        
        if (isEquil) {
            if (showAddRemove) {
                fprintf(pFileD, "\n\n************************************************");
                fprintf(pFileD, 
                    "\nADD RG %d TO EQUIL: Steps=%d t=%7.4e devious=%7.3e Rmin=%7.4e Rmax=%7.4e Ymin=%7.4e", 
                    RGn, totalTimeSteps, t, thisDevious, mineqcheck, maxeqcheck, Yminner);
            }
            
            for (int i = 0; i < niso; i++) {
                if (showAddRemove) {
                    fprintf(pFileD, "\n%s Z=%d N=%d Y=%8.4e Yeq=%8.4e Rprev=%8.4e Rnow=%8.5e",
                            isolabel[i], isoZ[i], isoN[i], isoY[i], isoYeq[i], eqcheck[i],
                            abs(isoY[i] - isoYeq[i]) / isoYeq[i]
                    );
                }
            }
            if (showAddRemove) 
                fprintf(pFileD, "\n************************************************\n");
            
        }
        
        // Set the activity array for each reaction in reaction group to true if not in 
        // equil and false if it is, if we are imposing equilibrium.
        
        if (doPE && t > equilibrateTime) {
            for (int i = 0; i < numberMemberReactions; i++) {
                int ck = memberReactions[i];
                reacIsActive[ck] = !isEquil;
            }
        }
            
    }    // End function computeEqRatios()
    
    
    // -----------------------------------------------------------
    // Function ReactionGroup::removeFromEquilibrium() to remove 
    // reaction group from equilibrium
    // -----------------------------------------------------------
    
    void removeFromEquilibrium() {
        
        isEquil = false;
        double thisDevious = abs((equilRatio - kratio) / kratio);
        
        if (showAddRemove) {
            fprintf(pFileD, "\n\n************************************************");
            fprintf(pFileD, 
                "\nREMOVE RG %d FROM EQUIL: Steps=%d t=%7.3e dt=%7.3e devious=%7.3e Rmin=%8.4e Rmax=%8.4e", 
                RGn, totalTimeSteps, t, dt, thisDevious, mineqcheck, maxeqcheck);
        }
        
        totalEquilRG --;
        
        for (int i = 0; i < niso; i++) {
            isEquil = false;
            if (showAddRemove) {
                fprintf(pFileD, "\n%s Z=%d N=%d Y=%8.4e Yeq=%8.4e Rprev=%8.4e Rnow=%8.5e",
                    isolabel[i], isoZ[i], isoN[i], isoY[i], isoYeq[i], eqcheck[i],
                    abs(isoY[i] - isoYeq[i]) / isoYeq[i]
                );
            }
        }
        
        for (int i = 0; i < numberMemberReactions; i++) {
            int ck = memberReactions[i];
            reacIsActive[ck] = true;         
            if (showAddRemove) {
                fprintf(pFileD, "\n Remove RG=%d %s RGflux=%7.4e flux[%d]=%7.4e", 
                    RGn, reaclabel[i], i, netflux, flux[i]);
            }
        }
        
        if(showAddRemove) 
            fprintf(pFileD, "\n************************************************\n");
    }
    
    
    // ----------------------------------------------------------------
    // Function ReactionGroup::speciesIsInRG() to determine if given 
    // species with index speciesIndex is in any of the reactions of a 
    // reaction group, where speciesIndex is the array index i for 
    // isotope quantitites like Z[i] or Y[i]. Returns true (1) if it
    // is and false (0) if not.
    // ----------------------------------------------------------------
    
    bool speciesIsInRG(int speciesIndex) { 
        
        // Loop over member reactions in the RG
        
        for(int i=0; i<numberMemberReactions; i++){

            // Loop over isotopes in reactants
            
            for (int j=0; j<numberReactants[i]; j++){ 
                if(isoindex[j] == speciesIndex){
                    return true;
                }
            }
            
            // Loop over isotopes in products
            
            for (int j=0; j<numberProducts[i]; j++){
                if(isoindex[j+numberReactants[i]] == speciesIndex){
                    return true;
                }
            }
        }
        
        return false;
        
    }       // End function speciesIsInRG(int)
    
};          // End class ReactionGroup



// ----------------------------------------------------------------
// Class Integrate with functions to integrate the reaction network.  
// Inherits from class Utilities.
// ----------------------------------------------------------------


class Integrate: public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions.  Its static functions can be called directly from the class
    // without having to instantiate.
    
    private:
        
    public:
        
        // Function to execute a single integration step.  Assumes that all fluxes have
        // already been calculated and relevant fluxes set to zero if PE approximation
        // and the reaction group has been judged to be in equilibrium.
        
        static double getTimestep(){

            // Store dt and sumX from previous step 
            dtLast = dt;
            sumXlast = sumX;
            
            //define variable for keeping track of any recalculations (mainly used for debugging purposes)
            int recountDown = 0;
            int recountUp = 0;
            
            //Define tolerances
            double const uptol = 1e-6;
            double const lowtol = 1e-10;
            double const tolC = 1e-8;
            double diffP;

            //define any dt adjusment parameters, dt. dtgrow, dtdec are global
            double dtnew;
            double dtmax;

            // start calculations that grow dt and update abundances + sumX
            dtnew = dtLast * dtgrow;
            dt = dtnew;
            dtmax = dt*1.10;

            updatePopulations(dt);
            sumX = Utilities::sumMassFractions();

            diffX = abs(sumX - 1.0);
            diffP = (abs(sumX - sumXlast)/sumX);

//Partial EQ addition
            if (doPE){
                if(mostDevious > deviousMax) {
                    dt = dt*dtdec;
                    recountDown++;

                    updatePopulations(dt);
                    sumX = Utilities::sumMassFractions();
                    diffX = abs(sumX - 1.0);

                     while(diffX > uptol){
                         dt = dt*dtdec;
                         recountDown++;

                         updatePopulations(dt);
                          sumX = Utilities::sumMassFractions();
                         diffX = abs(sumX - 1.0);
                        }
                }
                else if (mostDevious < deviousMin){

                    if (diffX < uptol){

                            if(diffP < tolC){
                                dt = dtmax;
                            }
                       else while (diffX < lowtol){
                        dt = dt*dtgrow;
                        recountUp++;

                        updatePopulations(dt);
                        sumX = Utilities::sumMassFractions();
                        diffX = abs(sumX - 1.0);

                        if(dt > dtmax){
                            dt = dtmax;
                            break;
                        }
                       }
                    }
                    else{
                        while(diffX > uptol){
                         dt = dt*dtdec;
                         recountDown++;

                         updatePopulations(dt);
                         sumX = Utilities::sumMassFractions();
                         diffX = abs(sumX - 1.0);
                        }
                    }
                }
            }

// without partial EQ
            else{

             while(diffX > uptol){
                dt = dt*dtdec;
                recountDown++;

                updatePopulations(dt);
                sumX = Utilities::sumMassFractions();
                diffX = abs(sumX - 1.0);

            }
            if (diffX < uptol){

               if(diffP < tolC){
                   dt = dtmax;
                }
               else while (diffX < lowtol){
                    dt = dt*dtgrow;
                    recountUp++;

                    updatePopulations(dt);
                    sumX = Utilities::sumMassFractions();
                    diffX = abs(sumX - 1.0);

                    if(dt > dtmax){
                        dt = dtmax;
                        break;
                    }
                }
            }
           }
        return dt;
        }


        // Function to do population update with current dt
        
        static void updatePopulations(double dtt){

            // If using the QSS approximation, apply QSS approximation to all isotopes
            
            if(doQSS){

                QSSupdate();
                
            }
            
            // If using asymptotic approximation, determine which species satisfy the
            // asymptotic condition.
            
            if(doASY){

                for(int i=0; i<ISOTOPES; i++){
                    isAsy[i] = checkAsy(keff[i]);
                }
                
                // Summarize results
                
                if(showAsyTest){
                    
                    bool asyck;
                    for(int i=0; i<ISOTOPES; i++){
                        asyck = false;
                        if(isAsy[i]) asyck=true;
                        double ck;
                        if( Y[i] > 0.0 ){
                            ck = FminusSum[i]*dtt/Y[i];
                        } else {
                            ck = 0.0;
                        }
                    }
                }
                
                // If Asy+PE, compute the matrix multiply for forward euler, with fluxes removed
                // by PE approximation (individual fluxes in RG that are in equilibrium) and 
                // asymptotic approximation (rows and columns of matrix)
                
                updateAsyEuler();
            }
            
        }
        
   
     
    // Function to update by the forward Euler method.  Returns the updated value of Y.
        
    static double eulerUpdate(int i, double fplusSum, double fminusSum, double y0, double dtt){
        
        double newY = y0 + (fplusSum-fminusSum)*dtt;
        
        if(diagnose2)
        fprintf(pFileD, 
        "\n  euler: %s t_i=%6.4e dt=%6.4e t_f=%6.4e k=%7.4e asycheck=%7.4e F+s=%6.4e F-=%6.4e dF=%6.4e Y0=%6.4e newY=%6.4e", 
        isoLabel[i], t, dtt, t+dtt, fminusSum/y0, fminusSum*dt/y0, fplusSum, fminusSum, fplusSum-fminusSum, y0, newY);

        return newY;     // New Y for forward Euler method
        
    }
    
    // Function to update by the asymptotic method using Sophia He formula. Returns
    // the updated value of Y.
    
    static double asymptoticUpdate(double fplus, double fminus, double y, double dtt){
        
        // Compute new Y for asymptotic method
        
        double newY = (y + fplus*dtt)/(1.0 + fminus*dtt/y);  
        
        if(diagnose2)
        fprintf(pFileD, 
        "\n  Asy: t_i=%6.4e dt=%6.4e t_f=%6.4e asycheck=%7.4e F+s=%6.4e F-=%6.4e dF=%6.4e Y0=%6.4e newY=%6.4e", 
        t, dtt, t+dtt, asycheck, fplus, fminus, fplus-fminus, y, newY);
        
        return newY;  
        
    }
    
    // Function Integrate::checkAsy() to determine whether an isotope satisfies the
    // asymptotic condition. Returns true (1) if it does and false (0) if not.
    // Two overloaded forms taking either Fminus and Y0, or keff as arguments.
    
    static bool checkAsy(double Fminus, double YY){
        
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
    
    
    /* Use the fluxes to update the populations for this timestep
     F or now we shall assume the asymptot*ic method. We determine whether each isotope 
     satisfies the asymptotic condition. If it does we update with the asymptotic formula. 
     If not, we update numerically using the forward Euler formula. */
    
    static void updateAsyEuler(){
        
        for(int i=0; i<numberSpecies; i++){ 
            
            if(isAsy[i]){
                Y[i] = asymptoticUpdate(FplusSum[i], FminusSum[i], Y0[i], dt);
            } else {
                Y[i] = eulerUpdate(i, FplusSum[i], FminusSum[i], Y0[i], dt);
            }
            
            X[i] = Y[i] * (double) AA[i];
            
        }
        
    }    // End function updateAsyEuler()
    
    
    // Function to update by the Quasi-Steady-State (QSS) approximation. This function 
    // replicates method steadyState(dt) in Java code. nitQSS is the number of
    // predictor-corrector iterations.  Found with Java code that increasing beyond 1
    // did not help much, so generally use nitQSS = 1, but code allows nitQSS > 1.
    
    static void QSSupdate(){
        
        // Iteration loop (nitQSS is number of predictor-corrector iterations; 
        // normally nitQSS=1)
        
        int nitQSS = 1;
        
        for (int i = 0; i < nitQSS; i++) {
            
            if(diagnoseQSS) fprintf(pFileD, "\n\nQSS_UPDATE t=%7.4e", t);
            
            ssPredictor();
            ssCorrector();
            
            // Recompute fluxes if another predictor-corrector iteration follows
            
            if (nitQSS > 1) {
                
                setReactionFluxes();
                
            }
        }
    }
    
    
    // ----------------------------------------------------------------------
    // Integrate::ssPredictor() to implement steady-state predictor step
    // ----------------------------------------------------------------------
    
    static void ssPredictor() {
        
        // Save current values of F+, F- and keff for later use. Y0[]
        // already contains the saved populations before update 
        // (i.e., from last timestep)
        
        for (int i = 0; i < ISOTOPES; i++) {

            FplusZero[i] = FplusSum[i];
            keffZero[i] = keff[i];
            
        }
        
        // Loop over all active isotopes and calculate the predictor
        // populations. Unlike for asymptotic method where we update with
        // the asymptotic approximation if kdt >= 1 and with forward euler if
        // if kdt < 1, we will update isotope abundances with the same 
        // QSS predictor-corrector, irrespective of value of keff*dt for 
        // an isotope.
        
        for (int i = 0; i < ISOTOPES; i++) {
            
                double kdt = keff[i] * dt;
                double deno = 1.0 + kdt*alphaValue(kdt);
                
                Y[i] = Y0[i] + (FplusSum[i] - FminusSum[i])*dt / deno;
                X[i] = Y[i] * (double)AA[i];
                
                if(diagnoseQSS) {
                    fprintf(pFileD, "\nSS PREDICTOR: %d t=%7.4e Fplus[i]=%7.4e",
                            i, t, FplusSum[i]);
                    fprintf(pFileD, " Fminus=%7.4e dt=%7.4e alph=%7.4e k=%7.4e",
                            FminusSum[i], dt, alphaValue(kdt), keff[i]);
                    fprintf(pFileD, "kdt=%7.4e deno=%7.4e Y=%7.4e",
                            kdt, deno, Y[i]);
                } 

        }
        
        // Update all fluxes for the corrector step
        
        setReactionFluxes();
        
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
                
                if(diagnoseQSS) {
                    
                    fprintf(pFileD, "\nSS CORRECTOR: %d t=%7.4e Fplus[i]=%7.4e",
                            i, t, FplusSum[i]);
                    fprintf(pFileD, " Fminus=%7.4e dt=%7.4e deno=%7.4e Y=%7.4e",
                            FminusSum[i], dt, 1 + alphaBar * kdt, Y[i]);
                    
                } 
                
                // For reference, keep track of the isotopes that would satisfy the
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


// Create an array of ReactionGroup objects RG[] to hold information and functions
// for the reactions groups in the network. Each element of the array will be
// a ReactionGroup object corresponding to a different reaction group of the network.
// Memory for array will be allocated dynamically below to the size
// given by numberRG (which is computed in the class ReactionVector)

ReactionGroup* RG;   // Dynamically allocated 1D array for reaction groups




// ------- Main CPU routine --------


int main() { 
    
    // Open a file for diagnostics output
    
    pFileD = fopen("gnu_out/diagnostics.data","w");
    
    // Write the time
    
    fprintf(pFileD, Utilities::showTime());
    printf("%s", Utilities::showTime());
    
    // Set labels and check consistency of choice for explicit algebraic methods set.
    // Generally we use either asymptotic (Asy) or quasi-steady-state (QSS) algorithms.
    // In either case we may choose to add the partial equilibrium (PE) algorithm. So
    // valid options are Asy, QSS, Asy+PE, and QSS+PE.
    
    if(doASY && !doPE){
        
       cout << "Using ASY method";
       fprintf(pFileD, "Using ASY method\n");
       doQSS = false;
       
    } else if (doQSS && !doPE) {
        
        cout << "Using QSS method";
        doASY = false;
        fprintf(pFileD, "Using QSS method\n");
        
    } else if (doASY && doPE){
        
        cout << "Using ASY+PE method";
        fprintf(pFileD, "Using ASY+PE method\n");
        
    } else if (doQSS && doPE){
        
        cout << "Using QSS+PE method";
        fprintf(pFileD, "Using QSS+PE method\n");
        
    }

    // Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
    // realistic calculation the temperature and density will be passed from the hydro 
    // code in an operator-split coupling of this network to hydro. Here we hardwire
    // them for testing purposes.  These will be used to calculate the reaction
    // rates in the network. If we assume operator splitting, the temperature 
    // and density are assumed constant for each network integration. But we allow
    // the possibility below to interpolate the temperature and density from a
    // hydrodynamical profile as a function of time.
    
    T9_start = 5.0;
    T9 = T9_start;
    rho_start = 1.0e8;
    rho = rho_start;
    
    // Initialize reacIsActive[] array to true;
    
    for (int i=0; i<SIZE; i++){ 
        
        reacIsActive[i] = true;
        
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
    
    // Print out some quantitites from the Reaction object reaction[].  
    
    for(int i=0; i<SIZE; i++){
        
        fprintf(pFileD, 
                "\n%d %s reacClass=%d reactants=%d products=%d isEC=%d isReverse=%d Q=%5.4f prefac=%5.4f", 
               reaction[i].getreacIndex(), 
               reaction[i].getreacChar(),  
               reaction[i].getreacClass(),
               reaction[i].getnumberReactants(),
               reaction[i].getnumberProducts(),
               reaction[i].getisEC(),
               reaction[i].getisReverse(),
               reaction[i].getQ(),
               reaction[i].getprefac()
        );
    }
    
    fprintf(pFileD, "\n\nReactantIndex[][] and ProductIndex[][]:\n\n");
    for(int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "%17s: ", reacLabel[i]);
        
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            fprintf(pFileD, "ReactantIndex[%d][%d]=%d ", i, j, ReactantIndex[i][j]);
        }
        
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            fprintf(pFileD, " ProductIndex[%d][%d]=%d", i, j, ProductIndex[i][j]);
        }
        
        fprintf(pFileD, "\n");
    }
    
    fprintf(pFileD, "\n\n\nREACLIB PARAMETERS FOR %d REACTIONS\n", SIZE);
    fprintf(pFileD, 
            "\n                                p0         p1         p2         p3         ");
    fprintf(pFileD, "p4         p5         p6");
    
    for (int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "\n%3d  %18s %10.4f", 
            i, reaction[i].getreacChar(), reaction[i].getp(0));
        
        for(int j=1; j<7; j++){
            fprintf(pFileD, " %10.4f", reaction[i].getp(j));
        }
        
    }
    
    fprintf(pFileD, "\n\nZ and N for reactants:\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "\n%3d %18s ", i, reaction[i].getreacChar());
        
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            fprintf(pFileD, " Z[%d]=%d", j, reaction[i].getreactantZ(j));
        }
        
        fprintf(pFileD, " ");
        
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            fprintf(pFileD, " N[%d]=%d", j, reaction[i].getreactantN(j));
        }
    }
    
    fprintf(pFileD, "\n\nZ and N for products:\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "\n%3d %18s ", i, reaction[i].getreacChar());
        
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            fprintf(pFileD, " Z[%d]=%d", j, reaction[i].getproductZ(j));
        }
        
        fprintf(pFileD, " ");
        
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            fprintf(pFileD, " N[%d]=%d", j, reaction[i].getproductN(j));
        }
    }

    fprintf(pFileD, 
    "\n\nreactantIndex for %d reactions (index of species vector for each reactant):\n", SIZE);

    for (int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "\n%d %18s ",i,reaction[i].getreacChar());
        
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            fprintf(pFileD, " reactantIndex[%d]=%d", j, reaction[i].getreactantIndex(j));
        }
    }
    
    fprintf(pFileD, 
    "\n\nproductIndex for %d reactions (index of species vector for each product):\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        
        fprintf(pFileD, "\n%d %18s ",i,reaction[i].getreacChar());
        
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            fprintf(pFileD, " productIndex[%d]=%d", j, reaction[i].getproductIndex(j));
        }
    }
    
    // Find the time intervals for plot output during the integration. After this
    // function is executed the plotSteps target time intervals for output will
    // be in the array plotTimeTargets[]. In the integration the ith output step will 
    // be triggered as soon as the time t is >= plottimeTargets[i].  The actual time of the output
    // (which will usually be slightly larger than plottimeTargets[i]) will be stored in tplot[i].
    // The variables start_time and stop_time define the range of integration.  The variable
    // startplot_time allows the plotting interval output in to be a subset of
    // the full integration interval.
    
    Utilities::log10Spacing(max(start_time, startplot_time), stop_time,
        plotSteps, plotTimeTargets);
    
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
    
    fprintf(pFileD, "\n\nNETWORK SPECIES VECTOR (%d components):\n\nIndex  Species    Z     N",
        numberSpecies);
    
    for(int i=0; i<numberSpecies; i++){
        fprintf(pFileD, "\n%5d    %5s  %3d  %4d", i, isoLabel[i], Z[i], N[i]);
    }
    
    fprintf(pFileD, "\n");
    
    // Use the information gleaned from ReactionVector::parseF() to define the reaction vectors
    // for the network using the static makeReactionVectors function of the class
    // ReactionVector.
    
    ReactionVector::makeReactionVectors();
    
    // Use static function ReactionVector::sortReactionGroups() to sort reactions into partial
    // equilibrium reaction groups by comparing reaction vectors.
    
    ReactionVector::sortReactionGroups();
    
    // Allocate dynamically an array of ReactionGroup objects of dimension numberRG, where
    // numberRG was determined by ReactionVector::sortReactionGroups() above.
    
    RG = (ReactionGroup*) malloc(sizeof(ReactionGroup)*numberRG);
    
    // Create ReactionGroup objects RG[] and assign values for various fields
    
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
    
    // Allocate 1D arrays to hold non-zero F+ and F- for all reactions for all isotopes,
    // the arrays holding the species factors FplusFac and FminusFac, and also arrays to hold 
    // their sums for each isotope. ReactionVector::parseF() must be run first because it determines 
    // totalFplus and totalFminus.
    
    Fplus = (double*) malloc(sizeof(double) * totalFplus);
    Fminus = (double*) malloc(sizeof(double) * totalFminus);
    FplusFac = (double*) malloc(sizeof(double) *totalFplus);
    FminusFac = (double*) malloc(sizeof(double) * totalFminus);
    FplusSum = (double*) malloc(sizeof(double) * numberSpecies);
    FminusSum = (double*) malloc(sizeof(double) * numberSpecies);
    
    // Allocate arrays that hold the index of the boundary between different isotopes in the
    // Fplus and Fminus 1D arrays. 
    
    FplusMax = (int*) malloc(sizeof(int) * numberSpecies);
    FplusMin = (int*) malloc(sizeof(int) * numberSpecies);
    FminusMax = (int*) malloc(sizeof(int) * numberSpecies);
    FminusMin = (int*) malloc(sizeof(int) * numberSpecies);

    // Allocate 1D arrays that will be used to map finite F+ and F- to the Flux array.
    
    FplusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
    FminusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
    FplusIsotopeIndex = (int*) malloc(sizeof(int) * totalFplus);
    FminusIsotopeIndex = (int*) malloc(sizeof(int) * totalFminus);
    
    // Allocate 1D arrays that will hold the index of the isotope for the F+ or F- term
    MapFplus = (int*) malloc(sizeof(int) * totalFplus);
    MapFminus = (int*) malloc(sizeof(int) * totalFminus);
    
    // Call function Reaction::setupFplusFminus() to set up F+ and F- index for each
    // isotope and to find non-vanishing F+ and F- source terms in network.
    // Function declared static so it can be called as Reaction::setupFplusFminus() 
    // without having to instantiate.
    
    Reaction::setupFplusFminus();
    
    
    
    // -----------------------------------------------
    // *** Begin main time integration while-loop ***
    // -----------------------------------------------
    
    printf("\n\n\n                 --- BEGIN TIME INTEGRATION ---\n");
    fprintf(pFileD, "\n\n\n\n                 --- BEGIN TIME INTEGRATION ---\n");
    
    dt = dt_start;              // Integration start time
    t = start_time;             // Current integration time
    totalTimeSteps = 1;         // Integration step counter
    totalEquilRG = 0;           // Number quilibrated reaction groups
    totalEquilReactions = 0;    // Number equilibrated reactions
    totalAsy = 0;               // Number asymptotic species
    ERelease = 0.0;             // Total E released from Q values
    int plotCounter = 1;        // Plot output counter
    fastestOverallRate = 0.0;   // Initialize fastest overall rate
    timeMaxRate = 0.0;          // Initialize slowest overall rate
    
    Utilities::startTimer();    // Start a timer for integration
    
    // Compute initial rates. If constant_T9 and constant_rho are true, rates won't
    // change in the integration and don't need to be computed again.  If either
    // T9 or rho change, the rates will be recomputed at each integration step.
    // Use functions of Reaction class to compute reaction rates. We have instantiated
    // a set of Reaction objects in the array reaction[i], one entry for each
    // reaction in the network. Loop over this array and call the computeRate()
    // function of Reaction on each object. 
    
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeConstantFacs(T9, rho);
        reaction[i].computeRate(T9, rho);
    }
    
    // Summarize computed rates
    
    fprintf(pFileD, "\nINITIAL COMPUTED RATES\n");
    for(int i=0; i<SIZE; i++){
        reaction[i].showRates();
    }
    
    if(constant_T9 && constant_rho){
        fprintf(pFileD, "\n\n**** Rates won't be computed again since T and rho won't ");
        fprintf(pFileD, "change in integration ****\n");
    } else {
        fprintf(pFileD, "\n\n**** Rates will be recomputed at each timestep since T and ");
        fprintf(pFileD, "rho may change ****\n");
    }
    

    
    while(t < stop_time){ 
        
        // Initialize fastest and slowest rates for this timestep
        
        fastestCurrentRate = 0.0;
        slowestCurrentRate = 1e30; 
        
        // Update Y0[] array with current Y[] values
        
        for(int i=0; i<ISOTOPES; i++){
            
            Y0[i] = Y[i];
            isotope[i].setY0(Y[i]);
            
        }
        
        // Specify temperature T9 and density rho. If constant_T9 = true, a constant
        // temperature is assumed for the entire network calculation, set by T9_start
        // above.  Otherwise (constant_T9 = false) we here interpolate the temperature
        // from a hydrodynamical profile for each timestep.  Likewise for the density.
        
        if(!constant_T9 && totalTimeSteps > 1){
            T9 = Utilities::interpolate_T(t);
        }
        if(!constant_rho && totalTimeSteps > 1){
            rho = Utilities::interpolate_rho(t);
        }
    
        // Use functions of Reaction class to compute reaction rates. We have instantiated
        // a set of Reaction objects in the array reaction[i], one entry for each
        // reaction in the network. Loop over this array and call the computeRate()
        // function of Reaction on each object. If constant_T9 and constant_rho are true, 
        // the rates only need be computed once as they won't change over this
        // integration. Otherwise they are recomputed for each integration step.
        
        if( (!constant_T9 || !constant_rho) && totalTimeSteps > 1){
            
            printf("**** RECOMPUTED RATES, timestep=%d\n",totalTimeSteps);
            
            for(int i=0; i<SIZE; i++){
                
                reaction[i].computeConstantFacs(T9, rho);
                reaction[i].computeRate(T9, rho);
                
            }
            
        }
        
        // Use functions of the Reaction class to compute fluxes.  We have instantiated
        // a set of Reaction objects in the array reaction[i], one entry for each
        // reaction in the network. Loop over this array and call the computeFlux()
        // function of Reaction on each object. Fluxes must be recomputed at each timestep
        // since they depend on the rates and the abundances. If temperature and density
        // are constant the rates won't change, but the fluxes generally will since
        // the abundances change even if the rates are constant as the network evolves.
        
        for(int i=0; i<SIZE; i++){
            
            reaction[i].computeFlux();
            
        }
        
        if(diagnose1){
            
            fprintf(pFileD, "\n\n\n--------- START NEW TIMESTEP: t_i = %7.4e Step=%d",
                    t, totalTimeSteps );
            fprintf(pFileD, " dt=%7.4e asyIsotopes=%d equilReaction=%d equilRG=%d ---------", 
                dt, totalAsy, totalEquilReactions, totalEquilRG );
            
        }
        
        // Call the static function Reaction::populateFplusFminus() to populate F+ and F-
        // for each isotope set up in setupFplusFminus() from master flux array computed
        // with setRGfluxes() above. We do it in this way because each flux is used in more than
        // one place but this way it is computed only once for each timestep.  This is efficient
        // because computing fluxes is the most time-consuming aspect of the algebraic explicit
        // algorithms.
        
        Reaction::populateFplusFminus();
        
        // Sum F+ and F- for each isotope
        
        Reaction::sumFplusFminus();
        
        // Summarize flux information if diagnose1=true
        
        if(diagnose1){
            
            fprintf(pFileD, "\n\nISOTOPE FLUXES:");
            
            for (int i=0; i<ISOTOPES; i++){
                
                fprintf(pFileD, 
                    "\n%d %s Y=%7.4e FplusSum=%7.4e FminusSum=%7.4e dF=%7.4e keff=%7.4e",
                    i, isotope[i].getLabel(), Y[i], //isotope[i].getY(), 
                    isotope[i].getfplus(), isotope[i].getfminus(), 
                    isotope[i].getfplus() - isotope[i].getfminus(),  
                    isotope[i].getkeff() );
                
            }
        }
        
        
        // Find max dY/dt and corresponding isotope for this timestep based on the initial
        // value of data for the timestep.
        
        getmaxdYdt();
        
        // Perform an integration step
        
        Integrate::getTimestep();
        
        // Update time to end of integration step and increment step counter.  The dt added here
        // has possibly been modified in doIntegrationStep() to satisfy tolerance conditions.
        
        t += dt; 
        totalTimeSteps ++; 


      if (totalTimeSteps > 100e6){
            Utilities::stopTimer();
            Utilities::plotOutput();
            return 0;
        }
        
        // Compute equilibrium conditions for the state at the end of this timestep (starting time
        // for next timestep) if partial equilibrium is being implemented.
        
        if(doPE && t > equilibrateTime){
            
            for(int i = 0; i < numberRG; i++) {
                RG[i].computeEquilibrium();
            }
            
            if(totalEquilRG > 0){
                
                if(diagnose2)
                fprintf(pFileD, 
                "\n\n********* BEGIN PE RESTORE: from t_i = %7.4e to t_f=%7.4e", t-dt, t);
                
                // Restore species in equilibrium to their unperturbed equilibrium values at 
                // the end of the timestep.  See the comments for function 
                // restoreEquilibriumProg() below for justification.
                
                restoreEquilibriumProg();
                
                if(diagnose2){
                    fprintf(pFileD, "\n\nISOTOPE FLUXES:");
                    for (int i=0; i<ISOTOPES; i++){
                        fprintf(pFileD, 
                            "\n%d %s Y=%7.4e F+Sum=%7.4e F-Sum=%7.4e dF=%7.4e keff=%7.4e",
                            i, isotope[i].getLabel(), Y[i], //isotope[i].getY(), 
                            isotope[i].getfplus(), isotope[i].getfminus(), 
                            isotope[i].getfplus() - isotope[i].getfminus(),  
                            isotope[i].getkeff());
                    }
                    
                    fprintf(pFileD, 
                    "\n\n********* END PE RESTORE: from t_i = %7.4e to t_f=%7.4e", t-dt, t);
                }
                
            }
        }
        
        // Count total asymptotic species
        
        totalAsy = 0;
        
        for(int i=0; i<ISOTOPES; i++){
            
            if (isAsy[i]){
                totalAsy ++;
            }
            
        }
        
        // Update the energy release variables based on Q values and fluxes for reactions
        
        double netdERelease = 0.0;
        
        for(int i=0; i<SIZE; i++){
            
            dERelease = reaction[i].getQ() * reaction[i].getflux();
            reaction[i].setdErate(dERelease);
            ERelease += (dERelease * dt);
            netdERelease += dERelease;
            
        }
        
        // ---------------------------------------------------------------------------------
        // Display and output to files updated quantities at plotSteps times corresponding
        // to (approximately) equally-spaced intervals in log_10(time). The target output 
        // times are generated by Utilities::log10Spacing() and stored in the array
        // plotTimeTargets[plotSteps]. A screen and plot file output is triggered if
        // t >= plotTimeTargets[plotCounter-1], so the actual output times may be slightly
        // larger than the target output times. Plots are made with respect to actual, not
        // target, output times.
        // ---------------------------------------------------------------------------------
        
        if(t >= plotTimeTargets[plotCounter-1]){
            
            // Optional output to diagnostic files
            
            if(showPlotSteps){
                fprintf(pFileD, "\n%s%s", dasher, dasher);
                fprintf(pFileD, 
                    "\n%d/%d steps=%d T9=%4.2f rho=%4.2e t=%8.4e dt=%8.4e asy=%d/%d sumX=%6.4f",
                    plotCounter, plotSteps, totalTimeSteps, T9, rho, t, dt, 
                    totalAsy, ISOTOPES, sumX);// diffX); // add diffX=%8.4e to print out diffX
                fprintf(pFileD, "\n%s%s", dasher, dasher);
                char tempest1[] = "\nIndex   Iso           Y           X        dY/dt";
                char tempest2[] = "        dX/dt           dY           dX\n";
                fprintf(pFileD, "%s%s", tempest1, tempest2);
                
                for(int i=0; i<ISOTOPES; i++){
                    fprintf(pFileD, "%5d %5s  %8.4e  %8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e\n", 
                        i, isoLabel[i], Y[i], X[i], isotope[i].getdYdt(), isotope[i].getdXdt(),
                        isotope[i].getdYdt()*dt, isotope[i].getdXdt()*dt
                    );
                }
                
                fprintf(pFileD, "%s%s\n", dasher, dasher);
            }
            
            // Output to plot arrays for this timestep
            
            tplot[plotCounter-1] = log10(t);
            dtplot[plotCounter-1] = log10(dt);
            
            // Log of E in erg and dE/dt in erg/g/s
            
            EReleasePlot[plotCounter-1] = log10( abs(ECON*ERelease) );
            dEReleasePlot[plotCounter-1] = log10( abs(ECON*netdERelease) );
            
            // Min and Max rate values (Rrate in sec^-1)
            
            slowestRatePlot[plotCounter-1] = slowestCurrentRate;
            slowestRateIndexPlot[plotCounter-1] = slowestCurrentRateIndex;
            fastestRatePlot[plotCounter-1] = fastestCurrentRate;
            fastestRateIndexPlot[plotCounter-1] = fastestCurrentRateIndex;
            
            sumXplot[plotCounter-1] = sumX;
            numAsyplot[plotCounter-1] = totalAsy;
            totalEquilRG = 0;
           // diffXplot[plotCounter -1] = log10(diffX); // add to output diffX to files for plotting
            
            for(int i=0; i<numberRG;i++){
                if(RG[i].getisEquil()) totalEquilRG ++;
            }
            
            numRG_PEplot[plotCounter-1] = totalEquilRG;
            
            for(int i=0; i<ISOTOPES; i++){
                
                Xplot[i][plotCounter-1] = X[i];
                FplusSumPlot[i][plotCounter-1] = isotope[i].getfplus();
                FminusSumPlot[i][plotCounter-1] = isotope[i].getfminus();
                
            }
            
            // Output to screen
            
            printf("\n%d/%d t=%7.4e dt=%7.4e Steps=%d Asy=%d/%d EqRG=%d/%d sumX=%5.3f diffX=%2.4e dE=%7.4e E=%7.4e",
                plotCounter, plotSteps, t, dt, totalTimeSteps, totalAsy, ISOTOPES, totalEquilRG, 
                numberRG, sumX, diffX, ECON*netdERelease, ECON*ERelease
            );//add diffX and diffX=%2.4e to plot diffX in output

            // Increment the plot output counter for next graphics output
            
            plotCounter ++;
        }
    
    }   // End time integration while-loop
    
    
    printf("\n\nEnd of integration");
    
    Utilities::stopTimer();      // Stop timer and print integration time

    // ------------------------------
    // *** End time integration ***
    // ------------------------------
    

    // Display abundances and mass fractions at end of integration

    printf("\nFINAL ABUNDANCES Y AND MASS FRACTIONS X\n");
    fprintf(pFileD, "\n\nFINAL ABUNDANCES Y AND MASS FRACTIONS X\n");

    for(int i=0; i<ISOTOPES; i++){
        
        printf("\n%d %s Y=%7.4e X=%7.4e", 
            i, isotope[i].getLabel(), Y[i], X[i]
        );
        
        fprintf(pFileD, "\n%d %s Y=%7.3e X=%7.3e", 
            i, isotope[i].getLabel(), Y[i], X[i]
        );
    }

    printf("\n\n");
    
    // Output of data to plot files after integration
    
    Utilities::plotOutput();
    
    
    // **************************************************
    // To test various function calls, insert code from 
    // functionTests.cpp here
    // **************************************************

    
    // Close output file
    
    fclose (pFileD);
   
    // Free allocated memory
    
    free(Fplus);
    free(Fminus);
    free(FplusFac);
    free(FminusFac);
    free(FplusSum);
    free(FminusSum);
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
// ************  FUNCTIONS USED IN MAIN  ********************
// **********************************************************


/* Function to adjust populations at end of timestep when in partial equilibrium 
 * to correct for deviations from equilibrium during the timestep. Uses progress 
 * variables. This function sets the progress variable for each reaction group to 
 * its equilibrium value at the end of the numerical timestep, and then sets the 
 * corresponding equilibrium values of other isotopes in the reaction group (since 
 * they are directly related to the equilibrium value of the progress variable for 
 * the reaction group). Then, the corrected abundances Y for all isotopes participating 
 * in partial equilibrium are averaged if they participate in more than one reaction 
 * group. Finally, after Ys are updated by their average equilibrium values, all 
 * abundances are rescaled so that total nucleon number is conserved by the overall 
 * timestep. Thus this function for restoring equilibrium does not require a matrix solution 
 * or Newton-Raphson iteration. It should scale approximately linearly with the network 
 * size. Contrast with the different approach taken in the method restoreEquilibrium()
 * of the Java code. */


void restoreEquilibriumProg() {

    int countConstraints = 0;
    int countEquilIsotopes = 0;
    int RGindy[numberRG];     // Array to hold index of RG in equilibrium
    sumX = Utilities::sumMassFractions();
    
    /* In general when we compute the equilibrium value of say alpha in the 
     *  reaction group alpha+16O <-> 20Ne, we are computing it using non-equilibrium 
     *  values of 16O and 20Ne (i.e., their values will not be the values that they 
     *  will have after this step). Add a while loop that permits iteration to try 
     *  to fix this. Preliminary tests indicate it has essentially no effect so set 
     * to one iteration for now. */
    
    int itcounter = 0;
    
    while (itcounter < 1) {
        
        countConstraints = 0;
        countEquilIsotopes = 0;
        itcounter ++;
        
        /* Compute equilibrium value of the Ys participating in equilibrium 
         *    starting from the value of Y at the end of the numerical timestep, 
         *    presently stored in Y[i]. Do so by first setting Y0[i] to the current 
         *    value of Y[i], which is the computed value at the END of the timestep. 
         *    Then evolve that initial value to the corresponding equilibrium value 
         *    algebraically by calculating the equilibrium value for that Y0[i] 
         *    and setting Y[i] to it (a form of operator splitting within the 
         *    network timestep). This work is done in evolveToEquilibrium(). */
        
        evolveToEquilibrium();
        
        // Inventory reaction groups in equilibrium
        
        for (int i = 0; i < numberRG; i++) {
            
            if ( RG[i].getisEquil() ) {
                
                RGindy[countConstraints] = i;
                countConstraints++;
                
                for (int j = 0; j < RG[i].getniso(); j++) {

                    int speciesIndy = RG[i].getisoindex(j);
                    
                    if (!isotopeInEquil[speciesIndy]) {
                        isotopeInEquil[speciesIndy] = true;
                        countEquilIsotopes++;
                    }
                }
            }
        }
        
        // Loop over reaction groups in equilibrium and compute equilibrated
        // Y[] averaged over all reaction groups that are in equilibrium and
        // contain the isotope.
        
        int numberCases;
        double Ysum;
        
        // Loop over isotopes checking for those in equilbrium in some RG
        
        for(int i=0; i<ISOTOPES; i++){
            
            // If isotope is in at least one equilibrated RG
            
            if (isotopeInEquil[i]) {
                numberCases = 0;
                Ysum = 0.0;
                
                // Find the equilibrated RGs the isotope appears in
                
                for(int j=0; j<numberRG; j++){
                    
                    if( RG[j].getisEquil() ){
                        
                        for(int k=0; k<RG[j].getniso(); k++){
                            if(i == RG[j].getisoindex(k)) {
                                Ysum += RG[j].getisoYeq(k);
                                numberCases ++;
                            }
                        }
                    }
                }
            }
            
            // Store Y for each isotope averaged over all reaction groups in 
            // which it participates
            
            Y[i] = Ysum/(double)numberCases;
            X[i] = Y[i]*(double)AA[i];
            
        }
    
    } // end while iteration loop
    
    
    // Set up renormalization of all Ys so that this integration step
    // conserves total particle number (sum of X = 1.0). Check mass 
    // fractions separately for isotopes participating in
    // equilibrium and those not.
    
    double sumXeq = Utilities::sumXEquil();
    double sumXNeq = Utilities::sumXNotEquil();
    
    // Factor to enforce particle number conservation
    
    double XcorrFac = 1.0 / (sumXeq + sumXNeq);

    // Loop over all isotopes and renormalize so sum X = 1
    
    for(int i=0; i<ISOTOPES; i++){
        
        X[i] *= XcorrFac;
        Y[i] = X[i] / (double)AA[i];
        Y0[i] = Y[i];
        
    }
    
    // Recompute total sumX 
    
    sumX = Utilities::sumMassFractions();
    
}


// ----------------------------------------------------------------------
// Function to set abundances for reaction groups in equilibrium to the
// current integrated value (end of timestep) and then evolve the 
// abundances algebraically to their equilibrium values.
// ----------------------------------------------------------------------

void evolveToEquilibrium() {

    for (int i = 0; i < numberRG; i++) {
        
        if ( RG[i].getisEquil() ) {
            
            // Loop over species in RG
            
            for (int j = 0; j < RG[i].getniso(); j++) {
                
                int indy = RG[i].getisoindex(j);
                Y0[indy] = Y[indy];
                RG[i].setisoY0(j, Y0[indy]);
                
            }
            
            // Compute equilibrium with new values of Y0
            
            RG[i].computeEquilibrium();
        }
    }
}


// ----------------------------------------------------------------------
// Function isoIsInRG(int isoindex, rgindex) to return 
// true if isotope labeled by species index isoindex is in the RG 
// labeled by rgindex and false otherwise. Not presently used
// since ReactionGroup::speciesIsInRG(speciesIndex) duplicates
// this functionality.
// ----------------------------------------------------------------------

bool isoIsInRG(int isoindex, int rgindex) {
    
    for (int j=0; j<RG[rgindex].getniso(); j++){
        
        if( RG[rgindex].getisoindex(j) == isoindex )
        return true;
        
    }
    
    return false;
}


// Find the maximum dY/dt for an isotope in the network

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



/* Function readNetwork to read the network data file line by line, with the filename as argument.
 * This file is expected to have 4 lines per isotope with the line structure
 * 
 *   isotopeSymbol A  Z  N  Y  MassExcess
 *   pf00 pf01 pf02 pf03 pf04 pf05 pf06 pf07
 *   pf10 pf11 pf12 pf13 pf14 pf15 pf16 pf17
 *   pf20 pf21 pf22 pf23 pf24 pf25 pf26 pf27
 * 
 * where isotopeSymbol is an isotope label, A=Z+N is the atomic mass number, Z is the proton number, 
 * N is the neutron number, Y is the initial abundance, MassExcess is the mass
 * excess in MeV, and the pf are 24 values of the partition function for that isotope at
 * different values of the temperature that will form a table for interpolation in temperature.
 * The assumed 24 values of the temperature for the partition function table are in array Tpf:
 * 
 * { 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 
 * 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 } in units of 10^9 K.
 * 
 * All fields on a line are separated by a blank space and there is no whitespace in the isotopeSymbol.
 * The type signature of these four lines corresponding to a single isotope is
 * 
 *  string int int int double double
 *  double double double double double double double double
 *  double double double double double double double double
 *  double double double double double double double double
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
 * A file with this format is written from the Java code to the file output/CUDAnetwork.inp using the
 * Java stream toCUDAnet.
 *  
 */

void readNetwork (char *fileName) {
    
    char line[60];
    char isoSymbol[5];
    int z, n, a;
    double y, mass;
    double pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7;
    
    // Open a file for reading 
    
    fr = fopen (fileName, "r");
    
    // Exit if the file doesn't exist or can't be read
    
    if( fr == NULL ){
        printf ("*** File Input Error: No readable file named %s\n", fileName);
        exit(1) ;
    }
    
    // Read in the file line by line
    
    int isoIndex = -1;
    int isoSubIndex = 3;
    
    if(displayInput==1) fprintf(pFileD, "\n--- Read in network and partition function ---\n");
    
    // Read lines until NULL encountered. Data lines can contain up to 60 characters.
    
    while(fgets(line, 60, fr) != NULL){
        
        isoSubIndex ++;
        
        if(isoSubIndex == 4){      // 1st line of data for first isotope
            
            isoSubIndex = 0;
            isoIndex ++;
            
            // Read 1st line
            
            sscanf (line, "%s %d %d %d %lf %lf", isoSymbol, &a, &z, &n, &y, &mass);
            
            if(displayInput == 1){
                fprintf(pFileD, "\n%s %d %d %d %f %f\n", isoSymbol, a, z, n, y, mass);
            }
            
            // Write data in 1st line to the Species object isotope[]
            
            isotope[isoIndex].setisoIndex(isoIndex);
            isotope[isoIndex].setisoLabel(isoSymbol);
            isotope[isoIndex].setZ(z);
            isotope[isoIndex].setN(n);
            isotope[isoIndex].setA(a);
            isotope[isoIndex].setY(y);
            isotope[isoIndex].setM(mass);
            
        } else {             // line contains partition function entries
            
            // Scan and parse a partition function line
            
            sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf", &pf0, &pf1, &pf2, &pf3, 
                &pf4, &pf5, &pf6, &pf7);
            
            if(displayInput == 1){
                fprintf(pFileD, "%f %f %f %f %f %f %f %f\n", 
                    pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7);
            }
            
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
        
        // Normally numberSpecies = ISOTOPES, but numberSpecies counts the 
        // actual number of reactions read in.
        
        numberSpecies = isoIndex + 1;  
        
    }
    
    // Close the file
    
    fclose(fr);
    
}    // End of function readNetwork (char *fileName)



/* Function to read rate parameter data file line by line, with filename as argument.
 * This file is expected to have one reaction per line with the line structure
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
    double p0, p1, p2, p3, p4, p5, p6, q, sf;
    int i0, i1, i2, i3, i4, i5, i6, i7;
    int ii[6];
    double tempp[7];
    int tempZN[4] = {-1, -1, -1, -1};
    
    // Open a file for reading  
    
    fr = fopen (fileName, "r");
    
    // Exit if the file doesn't exist or can't be read
    
    if( fr == NULL ){
        
        printf ("*** File Input Error: No readable file named %s\n", fileName);
        exit(1) ;
        
    }
    
    /* 
     * Read in the file line by line and parse into variables.  The expected
     * structure of each line is
     * 
     *     double double double double double double double string
     * 
     * each separated by a space, with no whitespace in the string.
     * (See http://stackoverflow.com/questions/2854488/reading-a-string-with-spaces-with-sscanf
     * for how to read string with spaces.)
     */
    
    int n = -1;
    int subindex = -1;
    
    if(displayInput == 1) fprintf(pFileD, "\n--- READ IN REACTIONS DATA---");
    
    // Read lines until NULL encountered. Lines can contain up to 120 characters.  In the
    // data file each reaction has 8 lines of entries.  The counter n holds the reaction number
    // and the counter subindex holds the line number for the current reaction.  The 
    // switch(subindex) determines which line we are reading (0-7) for a given reaction labeled by n.
    
    while(fgets(line, 120, fr) != NULL) {
        
        subindex ++;
        
        switch(subindex) {
            
            case 0:   // 1st line
                
                n++;
                
                sscanf (line, "%s %d %d %d %d %d %d %d %lf %lf %d", 
                    rlabel, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &sf, &q, &i7);
                
                // Store in the Reaction class instance reaction[n]
                
                ReactionPtr = &reaction[n];
                ReactionPtr->setreacIndex(n);
                ReactionPtr->setreacString(rlabel);   // setter also fills reacLabel[][] in main
                ReactionPtr->setreacGroupClass(i0);   // setter also fills RGclass[] in main
                ReactionPtr->setRGmemberIndex(i1);    // setter also fills RGMemberIndex[] in main
                ReactionPtr->setreacClass(i2);
                ReactionPtr->setnumberReactants(i3);  // setter also fills NumReactingSpecies[] in main
                ReactionPtr->setnumberProducts(i4);   // setter also fills NumProducts[] in main
                ReactionPtr->setisEC(i5);
                ReactionPtr->setisReverse(i6);
                ReactionPtr->setQ(q);
                ReactionPtr->setprefac(sf);
                ReactionPtr->setispeforward(i7);
                
                if(displayInput == 1){
                    
                    fprintf(pFileD, "\n\nReaction %d: ",n);
                    fprintf(pFileD, "%s reaclib=%d RG:(%s) RGclass=%d RGmemberIndex=%d",
                        reaction[n].getreacChar(), 
                        reaction[n].getreacIndex(),
                        reaction[n].getreacGroupChar(),
                        reaction[n].getreacGroupClass(), 
                        reaction[n].getRGmemberIndex());
                    fprintf(pFileD, 
                    "\n--------------------------------------------------------------------------------");
                    fprintf(pFileD, "\n%s %d %d %d %d %d %d %d %f %f", 
                        reacLabel[n], 
                        RGclass[n],
                        RGMemberIndex[n],
                        reaction[n].getreacClass(),
                        reaction[n].getnumberReactants(),
                        reaction[n].getnumberProducts(),
                        reaction[n].getisEC(),
                        reaction[n].getisReverse(),
                        reaction[n].getprefac(),
                        reaction[n].getQ()
                    );
                    
                }
                
                break;
                
            case 1:    // 2nd line
                
                sscanf (line, "%lf %lf %lf %lf %lf %lf %lf", &p0, &p1, &p2, &p3, &p4, &p5, &p6);
                
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
                
                if(displayInput == 1) fprintf(pFileD, "\n%f %f %f %f %f %f %f", 
                    reaction[n].getp(0),
                    reaction[n].getp(1),
                    reaction[n].getp(2),
                    reaction[n].getp(3),
                    reaction[n].getp(4),
                    reaction[n].getp(5),
                    reaction[n].getp(6)
                );
                
                break;
                
            case 2:    // 3rd line (Z of reactants)
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantZ(tempZN);    // setter also changes reacZ[][] in main
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        fprintf(pFileD, "\n  Reactant[%d]: Z=%d", mm, reacZ[n][mm]);
                    }
                    
                }
                
                break;
                
            case 3:   // 4th line (N of reactants)
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantN(tempZN);   // setter also changes reacN[][] in main
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        fprintf(pFileD, "\n  Reactant[%d]: N=%d", mm, reacZ[n][mm]);
                    }
                    
                }
                
                break;
                
            case 4:    // 5th line (Z of products)
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setproductZ(tempZN);    // setter also changes prodZ[][] in main
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        fprintf(pFileD, "\n  Product[%d]: Z=%d", mm, prodZ[n][mm]);
                    }
                    
                }
                
                break;
                
            case 5:    // 6th line (N of products)
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setproductN(tempZN);    // setter also changes prodN[][] in main
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        fprintf(pFileD, "\n  Product[%d]: N=%d", mm, prodN[n][mm]);
                    }
                    
                }
                
                break;
                
            case 6:    // 7th line (reactant species-vector index )
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberReactants(); k++){
                    tempZN[k] = ii[k];
                }
                
                ReactionPtr -> setreactantIndex(tempZN);   // setter also changes ReactantIndex[][]
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        fprintf(pFileD, "\n  ReactantIndex[%d]: N=%d", mm, ReactantIndex[n][mm]);
                    }
                    
                }
                
                break;
                
            case 7:    // 8th line (product species-vector index)
                
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                
                // Store in the Reaction class instance reaction[]
                
                ReactionPtr = &reaction[n];
                
                for(int i=0; i<4; i++){
                    tempZN[i] = -1;
                }
                
                for(int k=0; k<ReactionPtr -> getnumberProducts(); k++){
                    tempZN[k] = ii[k];
                }

                ReactionPtr -> setproductIndex(tempZN);    // setter also changes ProductIndex[][]
                
                if(displayInput == 1){
                    
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        fprintf(pFileD, "\n  ProductIndex[%d]: N=%d", mm, ProductIndex[n][mm]);
                    }
                    
                }
                
                subindex = -1;
                
                break;
                
        }

    }
    
    // Normally numberReactions=SIZE, but numberReactions counts the actual number of
    // reactions read in.
    
    numberReactions = n+1;
    
    fprintf(pFileD, "\n%d REACTIONS\n",numberReactions);
    
    // Close the file
    
    fclose(fr);           
    
}    // End of function readLibraryParams (char *fileName)



// Function to send the network isotopes, mass excesses, and the entries in the 
// partition function table for each isotope to the data file.

void writeNetwork() {

    fprintf(pFileD, "\n%d ISOTOPES IN NETWORK:\n\n",numberSpecies);
    fprintf(pFileD, "Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    
    for (int i=0; i<numberSpecies; i++){

        fprintf(pFileD, "%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n",  
            i, isoLabel[i], AA[i], Z[i], N[i], 
            Y[i], X[i], massExcess[i]);
        
    }
    
    // Write partition function table from isotope[] to data file
    
    fprintf(pFileD, "\n\nPARTITION FUNCTION TABLE from Species object isotope[]:\n");
    fprintf(pFileD, "\n T9 = ");
    
    for(int k=0; k<24; k++){
        fprintf(pFileD, "%4.2f ", isotope[0].getTpf(k));
    }
    
    for(int i=0; i<ISOTOPES; i++){

        fprintf(pFileD, "\n");
        fprintf(pFileD, "%-5s ",isotope[i].getLabel());
        
        for(int j=0; j<24; j++){ 
            fprintf(pFileD, "%4.2f ", isotope[i].getpf(j)); 
        }
        
    }

    fprintf(pFileD, "\n\n");
    
}   // End of function writeNetwork()



// Function to print out all the rates. The label can be used to distinguish cases
// if called more than once. NOT PRESENTRY USED.

void  writeRates(char *label) {
    
    printf("\n\nCOMPUTED RATES (%s):\n\n", label);
    
    for (int i=0; i<numberReactions; i++){
        
        printf("%d %s rate=%6.3e Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Q=%5.3f Prefac=%6.3e Reactants=%d\n",
            i, reaction[i].getreacChar(),  
            Rate[i]/reaction[i].getprefac(), Rate[i], 
            Y[reaction[i].getreactantIndex(0)],   
            Y[reaction[i].getreactantIndex(1)],
            Y[reaction[i].getreactantIndex(2)], 
            reaction[i].getQ(), 
            reaction[i].getprefac(), 
            reaction[i].getnumberReactants());
        
    }
    
}    // End of function writeRates(char *label)


// Helper function that can be called from another class to set some fields
// in the Reaction objects reaction[].

void setRG (int index, int RGclass, int RGindex) {

    reaction[index].setrgindex(RGindex);
    
    if(diagnose2)
    fprintf(pFileD, "\nsetRG: index=%d RGclass=%d RGindex=%d", 
        index, reaction[index].getreacGroupClass(), reaction[index].getrgindex());
}


// Helper function to create an array of ReactionGroup objects RG[i] and to set some fields 
// in the objects objects RG[].

void assignRG(){
    
    fprintf(pFileD, "\n\nREACTIONS IN RGclass[]:\n");
    
    for(int m=0; m<SIZE; m++){
        fprintf(pFileD, "\n%s RGclass[%d] = %d", reacLabel[m], m, RGclass[m]);
    }
    
    // Write out some fields for Reaction objects reaction[]
    
    fprintf(pFileD, "\n\n\nSOME FIELDS FOR THE %d Reaction OBJECTS reaction[]:\n", SIZE);
    
    for(int i=0; i<SIZE; i++){
        
        fprintf(pFileD, 
            "\nreaction[%d]: %s RGclass=%d #reac=%d #prod=%d RGmemberIndex=%d RG=%d",
            i, reaction[i].getreacChar(), 
            reaction[i].getreacGroupClass(),
            reaction[i].getnumberReactants(),
            reaction[i].getnumberProducts(),
            reaction[i].getRGmemberIndex(), 
            reaction[i].getrgindex()
        );
        
        int nummreac = reaction[i].getnumberReactants();
        int nummprod = reaction[i].getnumberProducts();
        
        // Write reactant symbols
        
        fprintf(pFileD, "\nRG=%d  REACTANTS: iso[0]=%s", 
        RG[i].getRGn(), isoLabel[reaction[i].getreactantIndex(0)]);
        
        if(nummreac > 1) fprintf(pFileD, " iso[1]=%s", isoLabel[reaction[i].getreactantIndex(1)]);
        if(nummreac > 2) fprintf(pFileD, " iso[2]=%s", isoLabel[reaction[i].getreactantIndex(2)]);
        
        // Write product Symbols
        
        fprintf(pFileD, 
        "  PRODUCTS: iso[%d]=%s", nummreac, isoLabel[reaction[i].getproductIndex(0)]);
        
        if(nummprod > 1) fprintf(pFileD, " iso[%d]=%s", nummreac+1, isoLabel[reaction[i].getproductIndex(1)]);
        if(nummprod > 2) fprintf(pFileD, " iso[%d]=%s", nummreac+2, isoLabel[reaction[i].getproductIndex(2)]);
        
        fprintf(pFileD, "\n");
        
    }
    
    // Loop to create and populate ReactionGroup objects RG[]
    
    fprintf(pFileD, "\n\nCREATING REACTION GROUPS RG[] AND POPULATING OBJECT FIELDS\n");
    
    for(int i=0; i<numberRG; i++){
        
        // Create array RG[i] of ReactionGroup objects
        
        RG[i] = ReactionGroup(i);
        RG[i].setnumberMemberReactions(RGnumberMembers[i]);
        
        // Set the reference reaction for the RG to be the first reaction in the RG
        // that has ifPEforward = true. This reference reaction will define the assumed
        // order of isotopes in the RG to be consistent with original Java code.
        
        int reffer = RG[i].setrefreac();

        int rgindex = -1;
        int upper1;
        int upper2;
        int ck1;
        
        // Loop over reactions of the network, picking out the members of RG[i] by the
        // condition that i = RGindex[j], where RGindex[j] holds the RG index for a
        // given reaction.
        
        for(int j=0; j<SIZE; j++){   
            
            // Following if() condition picks from the list of reactions one that is
            // in the ReactionGroup object RG[i]
            
            if(RGindex[j] == i){
                
                rgindex ++;          // Index for member reactions in RG
                
                RG[i].setmemberReactions(rgindex, j);
                RG[i].setrgclass(RGclass[j]);
                RG[i].setisForward(rgindex, isPEforward[j] );
                
                ck1 = RG[i].getmemberReactions(rgindex);  //reacIndex of member reaction in RG[i]
                
                RG[i].setnumberReactants(rgindex, reaction[ck1].getnumberReactants());
                RG[i].setnumberProducts(rgindex, reaction[ck1].getnumberProducts());
                RG[i].setreaclabel(rgindex, reaction[ck1].getreacString());
                
                upper1 = reaction[ck1].getnumberReactants();
                int rn = RG[i].getrefreac();
                int nrn = RG[i].getnumberReactants(rn);
                int ppp = RG[i].getmemberReactions(rn);
                
                // Loop over reactant isotopes within this reaction
                
                for(int k=0; k<nrn; k++){
                    
                    int qqq = reaction[ppp].getreactantIndex(k);
                    RG[i].setisoindex(k, qqq);
                    RG[i].setreactantIsoIndex(k, qqq);
                    RG[i].setisoZ(k, reaction[ppp].getreactantZ(k));
                    RG[i].setisoN(k, reaction[ppp].getreactantN(k));
                    RG[i].setisoA(k, reaction[ppp].getreactantA(k));
                    RG[i].setisolabel(k, isoLabel[qqq]);
                    
                }
                
                // Loop over product isotopes

                int nrn2 = RG[i].getnumberProducts(rn);
                upper2 = reaction[ck1].getnumberProducts();
                
                for(int k=0; k<nrn2; k++){
                    
                    int qqq = reaction[ppp].getproductIndex(k);
                    RG[i].setisoindex(k+nrn, qqq);
                    RG[i].setproductIsoIndex(k, qqq);
                    RG[i].setisoZ(k+nrn, reaction[ppp].getproductZ(k));
                    RG[i].setisoN(k+nrn, reaction[ppp].getproductN(k));
                    RG[i].setisoA(k+nrn, reaction[ppp].getproductA(k));
                    RG[i].setisolabel(k+nrn, isoLabel[qqq]);
                    
                }

                // Set the Ys in the RG

                int upk = RG[i].getnumberReactants(rn) + RG[i].getnumberProducts(rn);
                
                for(int k=0; k<upk; k++){
                    
                    int yindex = RG[i].getisoindex(k);
                    RG[i].setisoY(k, Y[yindex]);
                    
                }
                
                fprintf(pFileD, "\nreacIndex=%d memberIndex=%d %s RGclass=%d isForward=%d", 
                    RG[i].getmemberReactions(rgindex),
                    rgindex, RG[i].getreacString(rgindex),  
                    RG[i].getrgclass(),
                    RG[i].getisForward(rgindex)
                );
                
            }
        }
        
        // Populate the species index array for this RG using reference reaction reffer
        
        int RGclassRef = RGclass[RG[i].getmemberReactions(reffer)];
        RG[i].setniso(RGclassRef);
        fprintf(pFileD, "\nRG[%d]: refreac=%d RGclassRef=%d niso=%d Reactions=%d\n", 
               RG[i].getRGn(), RG[i].getrefreac(), RGclassRef, 
               RG[i].getniso(), RG[i].getnumberMemberReactions()
        );
    }
    
    // Summary of reaction groups
    
    for(int i=0; i<numberRG; i++){
        
        fprintf(pFileD, "\n\nSummary: RG=%d", RG[i].getRGn());
        int numr = RG[i].getnumberMemberReactions();
        
        for(int j=0; j<numr; j++){
            
            int reacID = RG[i].getmemberReactions(j);
            fprintf(pFileD, "\n%d %s iso[0]=%s iso[1]=%s iso[2]=%s iso[3]=%s", 
                j, reacLabel[reacID], RG[i].getisolabel(0),
                RG[i].getisolabel(1), RG[i].getisolabel(2),
                RG[i].getisolabel(3)
            );
            
        }
    }
    
    // Check that this function has assigned isotope indices in each
    // reaction group consistent with the order used for partial equilibrium
    // in the original Java code.
    
    fprintf(pFileD, "\n\n\nSUMMARY of order for isotopes for PE in ReactionGroup objects:");
    
    for(int i=0; i<numberRG; i++){
        
        int rn = RG[i].getrefreac();
        int upjj = RG[i].getnumberReactants(rn) + RG[i].getnumberProducts(rn);
        fprintf(pFileD, "\n\nRG=%d  RGclass=%d %s Species Index:", 
               i, RG[i].getrgclass(),
               Utilities::stringToChar( 
               reaction[RG[i].getmemberReactions(RG[i].getrefreac())].getreacGroupSymbol() )
        );
        
        for(int jj=0; jj<upjj; jj++){
            fprintf(pFileD, " iso[%d]=%d", jj, RG[i].getisoindex(jj));
        }
        
        fprintf(pFileD, "\n     ");
        for(int jj=0; jj<upjj; jj++){
            fprintf(pFileD, " isolabel[%d]=%s", jj, RG[i].getisolabel(jj));
        }
        
        fprintf(pFileD, "\n      REACTANTS: reactantIndex[0]=%d", RG[i].getreactantIsoIndex(0));
        for(int jj=1; jj<RG[i].getnumberReactants(rn); jj++){
            fprintf(pFileD, " reactantIndex[%d]=%d", jj, RG[i].getreactantIsoIndex(jj));
        }
        
        fprintf(pFileD, "\n      PRODUCTS: productIndex[0]=%d", RG[i].getproductIsoIndex(0));
        for(int jj=1; jj<RG[i].getnumberProducts(rn); jj++){
            fprintf(pFileD, " productIndex[%d]=%d", jj, RG[i].getproductIsoIndex(jj));
        }
        
        fprintf(pFileD, "\n      Z[0]=%d", RG[i].getisoZ(0));
        for(int jj=1; jj<upjj; jj++){
            fprintf(pFileD, " Z[%d]=%d", jj, RG[i].getisoZ(jj));
        }
        
        fprintf(pFileD, "\n      N[0]=%d", RG[i].getisoN(0));
        for(int jj=1; jj<upjj; jj++){
            fprintf(pFileD, " N[%d]=%d", jj, RG[i].getisoN(jj));
        }
        
        fprintf(pFileD, "\n      A[0]=%d", RG[i].getisoA(0));
        for(int jj=1; jj<upjj; jj++){
            fprintf(pFileD, " A[%d]=%d", jj, RG[i].getisoA(jj));
        }
        
        fprintf(pFileD, "\n      isoY[0]=%8.5e", RG[i].getisoY(0));
        for(int jj=1; jj<upjj; jj++){
            fprintf(pFileD, " isoY[%d]=%8.5e", jj, RG[i].getisoY(jj));
        }
        
        fprintf(pFileD, "\n      isoYeq[0]=%8.5e", RG[i].getisoYeq(0));
        for(int jj=1; jj<upjj; jj++){
            fprintf(pFileD, " isoYeq[%d]=%8.5e", jj, RG[i].getisoYeq(jj));
        }
    }
    
}       // End function assignRG()


// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesfplus(int index, double fp){
    
    isotope[index].setfplus(fp);
    
}

// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesfminus(int index, double fm){
    
    isotope[index].setfminus(fm);
    isotope[index].setkeff(fm/(isotope[index].getY0()+1.0e-30));
    keff[index] = isotope[index].getkeff();
    
}

// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesdYdt(int index, double dydt){
    
    isotope[index].setdYdt(dydt);
    
}

// Helper function that can be called from another class to set
// flux fields in the reaction[] objects

void setReactionFluxes(){
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeFlux();
    }
    
    Reaction::populateFplusFminus();
    Reaction::sumFplusFminus();
}



