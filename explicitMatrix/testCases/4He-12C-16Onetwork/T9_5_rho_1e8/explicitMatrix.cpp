/*
 * Code to implement explicit algebraic integration of astrophysical thermonuclear networks.
 * Execution assuming use of Fedora Linux: Compile with
 * 
 *     gcc explicitMatrix.cpp -o explicitMatrix -lgsl -lgslcblas -lm -lstdc++
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
 * AUTHORS:
 * ---------------
 * Nick Brey
 * Ashton DeRousse
 * Adam Cole
 * Amelia Konomos
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

#define ISOTOPES 3                    // Max isotopes in network (e.g. 16 for alpha network)
#define SIZE 8                        // Max number of reactions (e.g. 48 for alpha network)
#define plotSteps 10                  // Number of plot output steps

#define LABELSIZE 35                  // Max size of reaction string a+b>c in characters
#define PF 24                         // Number entries in partition function table for each isotope
#define THIRD 0.333333333333333
#define TWOTHIRD 0.66666666666667

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
#define PRINT_CPU (printf("Timer: %g ms used", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define PRINT_CPU_TEST (printf("\nTimer Test: %g ms used by CPU\n", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));

// File pointer for data read-in
FILE *fr;            

// Filename for input rates library data. The file rateLibrary.data output by the Java code through 
// the stream toRateData has the expected format for this file.  Standard test cases: 
// rateLibrary_alpha.data, rateLibrary_150.data, rateLibrary_365.data, rateLibrary_nova134.data,
// rateLibrary_3alpha.data, rateLibrary_pp.data.

char rateLibraryFile[] = "data/rateLibrary_3alpha.data";  

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format for 
// this file. Standard test cases: CUDAnet_alphasolar.inp, CUDAnet_150solar.inp,
// CUDAnet_365solar.inp, CUDAnet_nova134.inp, CUDAnet_3alpha.inp, CUDAnet_pp.inp.

char networkFile[] = "data/CUDAnet_3alpha.inp";

// Control diagnostic printout of details (1 to print, 0 to suppress)
static const int displayInput = 0;
static const int showParsing = 0;
static const int showFparsing = 0;
static const int showFluxCalc = 0; //1;
static const int showRVdetails = 0;
static const int showRGsorting = 0;
static const int showAsyTest = 0; //1;
static const int showFunctionTests = 0;
// Whether to write message when RG added/removed from equil
bool showAddRemove = true; 

// Function Signatures:
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

// Control which explicit algebraic approximations are used. Eventually
// this should be set from a data file. To use asymptotic set doASY true
// (which toggles doQSS to false). To use quasi-steady-state (QSS), set 
// doASY false (which toggles doQSS to true). doPE can be true or false 
// with either Asymptotic or QSS.

bool doASY = true;            // Whether to use asymptotic approximation
bool doQSS = !doASY;          // Whether to use QSS approximation 
bool doPE = false;             // Implement partial equilibium also

// String holding integration method in use.  Possibilities are
// ASY, QSS, ASY+PE, QSS+PE

string methstring; 

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

// Array to hold whether given species satisfies asymptotic condition
// True (1) if asyptotic; else false (0).

bool isAsy[ISOTOPES];

// Force a constant timestep constant_dt for testing purposes by
// setting constantTimestep=true.  Normally constantTimestep=false
// for adaptive timestepping.

bool constantTimestep = false;    // Adaptible timestep if false
double constant_dt = 1.1e-9;     // Value of constant timestep

// Integration time data.  Start and stop times hardwired for testing
// here but in applications they would be variables supplied by
// calling programs. Likewise for dt_start.

double start_time = 1.0e-12;         // Start time for integration
double logStart = log10(start_time); // Base 10 log start time
double stop_time = 1.0e-11;           // Stop time for integration
double logStop = log10(stop_time);   // Base-10 log stop time
double dt_start = 1.0e-13;           // Initial value of integration dt
double dt;                           // Current integration timestep
double t;                            // Current time in integration
int totalTimeSteps;                  // Number of integration timesteps taken
double deltaTime;                    // dt for current integration step
int totalAsy;                        // Total number of asymptotic isotopes

double sumX;                         // Sum of mass fractions X(i).  Should be 1.0.
double diffX;                        // sumX - 1.0

double Rate[SIZE];       // Computed rate for each reaction from Reactions::computeRate()
double Flux[SIZE];       // Computed flux for each reaction from Reactions::computeFlux()


// --- Species data in the following arrays also contained in fields of class Species

int Z[ISOTOPES];                 // Array holding Z values for isotopes
int N[ISOTOPES];                 // Array holding N values for isotopes
int AA[ISOTOPES];                // Array holding A values for isotopes
double Y[ISOTOPES];              // Array holding abundances Y for isotopes
double X[ISOTOPES];              // Array holding mass fractions X for isotopes
double massExcess[ISOTOPES];     // Array holding mass excesses for isotopes
char isoLabel[ISOTOPES][5];      // Isotope labels (max length 5 characters; e.g. 238pu)

// -----------


// --- Reaction data in the following arrays also contained in fields of class Reaction

char reacLabel[SIZE][LABELSIZE]; // Char array of reaction labels (e.g. he4+c12-->o16) 
int RGclass[SIZE];               // Reaction Group class (PE) for reaction (1-5)
int RGMemberIndex[SIZE];         // Member index within its reaction group
string RGstring[SIZE];           // Schematic RG; e.g. a <->b+c
int NumReactingSpecies[SIZE];    // Number of reactant isotopes for the reaction
int NumProducts[SIZE];           // Number of product isotopes for the reaction
int reacZ[SIZE][4];              // Holds Z for each reactant isotope
int reacN[SIZE][4];              // Holds N for each reactant isotope
int prodZ[SIZE][4];              // Holds Z for each product isotope
int prodN[SIZE][4];              // Holds N for each product isotope
int ReactantIndex[SIZE][4];      // Index of isotope vector for each reactant isotope
int ProductIndex[SIZE][4];       // Index of isotope vector for each product isotope
int isPEforward[SIZE];           // Whether labeled "forward" reaction in PE scheme

// -----------


int numberSpecies;                // Actual # species in network (generally = ISOTOPES)
int numberReactions;              // Actual # reactions in network (generally = SIZE)

// Array with entries +1 if a reaction increases the population of the isotope (contributes to 
// F+), -1 if it decreases it (contributes to F-) and 0 if the reaction does not change the population
// of the isotope. This array is populated in the function ReactionVector::parseF().  It is 
// characteristic of the structure of the network and thus has to be calculated only once for a 
// given network.

int reacMask[ISOTOPES][SIZE]; 

// Define an array rv[] and corresponding pointers that will hold GSL vectors corresponding
// to the reaction vectors for the system.  This will be implemented in the function
// makeReactionVectors() of the class ReactionVector.

gsl_vector rv[SIZE];   // Array of type gsl_vector to hold GSL vectors
gsl_vector *rvPt;      // Pointer to rv[] array

// Total number of F+ and F- terms in the network
int totalFplus = 0;
int totalFminus = 0;

// Arrays to hold non-zero fluxes in the network. Corresponding memory will be allocated dynamically 
// below with malloc

double* Fplus;        // Dynamically-allocated 1D array for non-zero F+ (Dim totalFplus)
double* Fminus;       // Dynamically-allocated 1D array for non-zero F- (Dim totalFminus)

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

string tempest;   // Utility string to hold temporary quantities
string dasher = "---------------------------------------------";


// ---------------------------------
// Partial equilibrium quantities
// ---------------------------------

int numberRG;          // Number of partial equilibrium reaction groups

// Define array to hold the reaction group index for each reaction. There are n reaction
// groups in the network and each reaction belongs to one reaction group.  RGindex[m] is
// the index (0, 1, ... n) of the reaction group to which reaction m belongs. 
// This array is populated by the function ReactionGroup::sortReactionGroups()

int RGindex[SIZE];

double mostDevious = 0.0;     // Largest deviation of equilibrium k ratio from equil
int mostDeviousIndex;         // Index of RG with mostDevious
double maxDevious = 0.5;      // Max allowed deviation of Y from equil value in timestep

// Whether to compute and display partial equilibrium quantities
bool equilibrate = true;

// Whether actually to impose partial equilibrium
bool imposeEquil = true;  

// Time to begin trying to impose partial equilibrium.  Hardwired for now, but eventually
// this should be determined by the program.  In the Java version this was sometimes
// needed because starting PE test too early could lead to bad results.  This is 
// probably a coding error in the Java version, since if operating properly nothing should
// be changed at a timestep if nothing satisfies PE condition.  Thus, we should not need
// this in a final version for stability, but it might still be useful since early in
// a calculation typically nothing satisfies PE, so checking for it is a waste of time.
// However, check should not be costly.

double equilibrateTime = 1.0e-9; 

double equiTol = 0.01;      // Tolerance for checking whether Ys in RG in equil 

double Yminner;             // Current minimum Y in reaction group
double mineqcheck;          // Current minimum value of eqcheck in reaction group
double maxeqcheck;          // Current max value of eqcheck in reaction group

bool reacIsActive[SIZE];    // False if reaction has been removed by PE

int totalEquilReactions;    // Total equilibrated reactions for isotope
int totalEquilRG;           // Total equilibrated reaction groups

// Threshold abundance for imposing equil in reactions.  There may be numerical
// issues if the PE algorithm is imposed for very small abundances early in
// the calculation.

double Ythresh = 0.0;       

gsl_matrix *fluxes;
gsl_vector *abundances;


// -------------------------------------------------------------------------
// Arrays to hold output quantities at plot timesteps. The variable
// plotSteps is hardwired above, but eventually should be input
// at time of execution.
// -------------------------------------------------------------------------

double plotTimeTargets[plotSteps];     // Target plot times for plot step

double tplot[plotSteps];               // Actual time for plot step
double dtplot[plotSteps];              // dt for plot step
double Xplot[ISOTOPES][plotSteps];     // Mass fractions X
double sumXplot[plotSteps];            // Sum of mass fractions
double numAsyplot[plotSteps];          // Number asymptotic species
double numRG_PEplot[plotSteps];        // Number RG in PE
double EReleasePlot[plotSteps];        // Integrated energy release
double dEReleasePlot[plotSteps];       // Differential energy release

// Following control which mass fractions are exported to the plotting
// file.  The entries in plotXlist[] are the species indices for the
// isotopes in the network to be plotted.

int plotXlist[] = {1, 2, 3};           // Species index for X to be plotted 
int LX;                                // Length of plotXlist array




//----------------CLASS DEFINITIONS ----------------


/* Class Utilities to hold utility useful utility functions.  Functions are
 * declared static so that they can be invoked without having to instantiate
 * objects of type Utilities.  For example, Utilities::returnNetIndexZN (Z, N).
 * We will often have other classes inherit from Utilities so that they can
 * access its static methods.  This class definition must precede definitions of 
 * classes that inherit from it.
 */


class Utilities{
    
    private:
    
    public:
        
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
            
            // Open a file for ascii output. Assumes that the subdirectory
            // gnu_out already exists. If it doesn't, will compile but
            // crash when executed.
            
            FILE * pFile;
            pFile = fopen("gnu_out/gnufile.data","w");
            
            // Get length of array plotXlist holding the species indices for isotopes
            // that we will plot mass fraction X for.
            
            LX = sizeof(plotXlist)/sizeof(plotXlist[0]);
            
            string str1 = "#    t     dt   |E| |dE/dt|  Asy  Equil  sumX";
            string app = "  ";
            string app1;
            string Xstring = "X(";
            string iso;
            
            fprintf(pFile, "# %s:  %d integration steps\n",
                stringToChar(methstring), totalTimeSteps
            );
            
            fprintf(pFile, "# All quantities except Asy, RG_PE, and sumX are log10(x)\n");
            fprintf(pFile, "# Log of absolute values for E and dE/dt as they can be negative\n");
            fprintf(pFile, "# Units: t and dt in s; E in erg; dE/dt in erg/g/s; others dimensionless \n");
            fprintf(pFile, "#\n");
            
            // Write header for gnuplot file
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[i];
                app.append(Xstring);
                app.append(iso);
                app.append(")     ");
            }
            str1.append(app);
            str1.append("\n");
            fprintf(pFile, stringToChar(str1));
            printf("\n");
            
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
                );
                
                // Now add one data field for each X(i) in plotXlist[]. Add
                // 1e-24 to X in case it is identically zero since we are
                // taking the log.
                
                for(int j=0; j<LX; j++){
                    fprintf(pFile, " %5.3e", log10(Xplot[j][i]+1e-24));
                }
                
                fprintf(pFile, "\n");
            }
            
            fclose (pFile);
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::returnSumX() to return the current sum of the
        // mass fractions X(i) in the network. If the network conserves particle
        // number this sum should be equal to 1.0.
        // -------------------------------------------------------------------------
        
        static double returnSumX(void) {
            double sum = zerod;
            for(int i=0; i<ISOTOPES; i++){
                sum += X[i];
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
        // to printf typically displays garbage because of type issues in the 
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
        double YY;           // abundance  Y = X/A
        double XX;           // mass fraction  X = Y*A
        double MassExcess;   // mass excess
        double pf[24];       // partition function entries
        double fplus;        // Total flux currently adding to abundance of this isotope
        double fminus;       // Total flux currently adding to abundance of this isotope
        double dYdt;         // Current dY/dt for this isotope
        double dXdt;         // Current dX/dt for this isotope
        
        // Temperatures in units of 10^9 K for partition function table (see pf[]). 
        const double Tpf[24] = { 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 
            0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 
            4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f };
    
    public:
        
        // Public setter functions to set values of private class fields
        
        void setisoIndex(int iso){
            isoindex = iso;
        }
        
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
        
        void setdYdt(double d){
            dYdt = d; 
            dXdt = d*(double)A;
        }
        
        void setdXdt(double d){
            dXdt = d;
            dYdt = d/(double)A;
        }
        
        
        // Public getter functions to return values of class fields
        
        int getZ() {return ZZ; };
        int getN() {return NN; };
        int getA() {return A; };
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
        
        double getfminus(){
            return fminus;      
        }
        
        double getfplus(){
            return fplus; 
        }
        
        double getdYdt(){
            return dYdt;              
        }
        
        double getdXdt(){
            return dXdt;              
        }
        
        
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
        //      strcpy(cs, &s[0]);	// or strcpy(cs, s.c_str());
        //      printf("\n\nstring=%s\n", strcpy(cs, &s[0]));
        //
        // The function getreacChar() below returns the string reacString as a
        // pointer to a character array that will work in printf. Alternatively,
        // Utilities::stringToChar() will do same thing.
        
        string reacGroupClassLett;   // Letter equivalent (A-E) for reacGroupClass
        string reacGroupSymbol;      // Schematic equil reaction (e.g. a+b<->c)
        int numberReactants;         // Number species on the left side of reaction
        int numberProducts;          // Number species on the right side of reaction
        string reacString;           // String describing reaction
        string resonanceType;        // Whether resonant (r) or non-resonant (nr)
        int isEC;                    // Whether electron capture reaction (1) or not (0)
        int isReverse;               // Whether reverse reaction (1) or not (0)
        int ispeforward;             // Whether reactions is labeled "forward" in PE scheme
        
        double Q;                    // Q-value for reaction
        double prefac;               // The eta prefac for rates
        double p[7];                 // ReacLib parameters
        int reactantZ[3];            // Array holding Z of reactants
        int reactantN[3];            // Array holding N of reactants
        int productZ[4];             // Array holding Z of products
        int productN[4];             // Array holding N of products
        int reactantIndex[3];        // Index of species isotope vector for each reactant isotope
        int productIndex[4];         // Index of species isotope vector for each product isotope
        
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
        double flux;                 // Current flux of reaction
        char cs[20];                 // Utility character array
        char ccs[20];                // Utility character array
        
        string ss;                   // Utility string
  
  
    public:
        
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
        
        void setdensfac(double d){ densfac = d;}
        
        void setrate(double r){ rate = r; }
        
        void setRrate(double r){ Rrate = r; }
        
        void setflux(double f){ flux = f; }
        
        
        // Public Reaction getter methods to get values in private fields
        
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
        
        double getrate(){ return rate; }
        
        double getRrate(){ return Rrate; }
        
        double getflux(){ return flux; }
        
        
        // Function Reaction::setupFplusFminus() to set up F+ and F- index for each
        // isotope and to find non-vanishing F+ and F- source terms in network.
        // Method declared static so it can be called as Reaction::setupFplusFminus() 
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
                printf("\n\n--- MAX F+ and F- INDEX FOR EACH ISOTOPE ---\n");	
                for(int i=0; i<numberSpecies; i++)
                {
                    printf("\nIsotope index = %d  %s  Max index F+ = %d  Max index F- = %d", 
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
                printf("\n\n\n--- %d NON-VANISHING F+ SOURCE TERMS ---\n", totalFplus);
                printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
                       isoLabel[FplusIsotopeIndex[0]], FplusIsotopeIndex[0],
                       numFluxPlus[FplusIsotopeIndex[0]]);
                for(int i=0; i<totalFplus; i++){
                    printf("\n   Isotope index = %d F+ index = %d Reac index = %d  %s", 
                           FplusIsotopeIndex[i], i, MapFplus[i], 
                           reacLabel[MapFplus[i]]); 
                    if(i == (FplusIsotopeCut[FplusIsotopeIndex[i]] - 1)  && i != totalFplus-1)
                    {
                        printf("\n");
                        printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
                               isoLabel[FplusIsotopeIndex[i+1]], FplusIsotopeIndex[i+1],
                               numFluxPlus[FplusIsotopeIndex[i+1]]);
                    }
                }	
                printf("\n\n\n--- %d NON-VANISHING F- SOURCE TERMS ---\n", totalFminus);
                printf("\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):", 
                       isoLabel[FminusIsotopeIndex[0]], FminusIsotopeIndex[0],
                       numFluxMinus[FminusIsotopeIndex[0]] );
                for(int i=0; i<totalFminus; i++){
                    printf("\n   Isotope index = %d F- index = %d Reac index=%d  %s", 
                           FminusIsotopeIndex[i], i, MapFminus[i], reacLabel[MapFminus[i]]);
                    if(i == (FminusIsotopeCut[FminusIsotopeIndex[i]] - 1) && i != totalFminus-1 ){
                        printf("\n");
                        printf("\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):", 
                               isoLabel[FminusIsotopeIndex[i+1]], FminusIsotopeIndex[i+1],
                               numFluxMinus[FminusIsotopeIndex[i+1]]
                        );
                    }
                }
                printf("\n\n");
            }
        }
        
        
        // Function Reaction::populateFplusFminus() to populate F+ and F- for each
        // isotope set up in setupFplusFminus() from master flux array. Method declared
        // static so it can be called as Reaction::populateFplusFminus() without having
        // to instantiate.
        
        static void populateFplusFminus(){
            
            // Populate the F+ and F- arrays from the master Flux array
            
            //printf("\n\nVALUES F+\n\n");

            for(int i=0; i<totalFplus; i++){
                int indy = MapFplus[i];
                Fplus[i] = FplusFac[i]*Flux[indy];
//                 printf("i=%d FplusFac=%3.1f Flux=%7.3e Fplus=%7.3e\n", 
//                     i, FplusFac[i], Flux[indy], Fplus[i]);
            }
            
            //printf("\nVALUES F-\n\n");
            
            for(int i=0; i<totalFminus; i++){
                int indy = MapFminus[i];
                Fminus[i] = FminusFac[i]*Flux[indy];
//                 printf("i=%d FminusFac=%3.1f Flux=%7.3e Fminus=%7.3e\n", 
//                        i, FminusFac[i], Flux[indy], Fminus[i]);
            }
            //printf("\n");
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
        
        // Reaction::computeRate(double, double) to compute rates at T and rho. The quantity rate
        // is the temperature-dependent part.  The quantity Rrrate is rate multiplied by
        // appropriate density and statistical factors, which give units of s^-1.  The 
        // flux follows from multiplying Rrate by appropriate abundances Y in computeFlux().
        
        void computeRate(double T9, double rho){
            
            // Temperature-dependent rate from ReacLib library
            rate = expf( p[0] + t1*p[1] + t2*p[2] + t3*p[3] + t4*p[4] + t5*p[5] + t6*p[6] );
            setrate(rate);

            // Full rate factor in s^-1 (rate above multiplied by density factors)
            Rrate = getdensfac() * rate;

            setRrate(Rrate);
            
            // Write to rate array in main
            Rate[getreacIndex()] = Rrate;
            
            printf("\n%d %19s densfac=%6.3e rate= %6.3e Rrate=%6.3e", 
                   getreacIndex(), getreacChar(), getdensfac(), getrate(), getRrate());
            
        }
        
        
        // Reaction::computeFlux() to compute the current flux for reaction corresponding 
        // to this Reaction object.
        
        void computeFlux(){
            
            string s;
            
            switch(numberReactants){
                
                case 1:    // 1-body reactions
                    
                    flux = Rrate*Y[ reactantIndex[0] ];	
                    Flux[getreacIndex()] = flux;         // Put in flux array in main
                    if(showFluxCalc == 1){
                        printf("\n%d %18s reactants=%d iso0=%d Rrate=%7.3e Y1=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0],  Rrate, 
                            Y[ reactantIndex[0] ], flux);
                    }
                    break;
                    
                case 2:	   // 2-body reactions	
                    
                    flux = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ]; 	
                    Flux[getreacIndex()] = flux;         // Put in flux array in main
                    if(showFluxCalc == 1){
                        s = "\n%d %18s reactants=%d iso0=%d iso1=%d Rrate=%7.3e ";
                        s += "Y1=%7.3e Y2=%7.3e Flux=%7.3e";
                        printf(Utilities::stringToChar(s),
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], 
                            reactantIndex[1], Rrate, Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], flux);
                    }
                    break;
                    
                case 3:	   // 3-body reactions
                    
                    flux = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ] * Y[ reactantIndex[2] ];
                    Flux[getreacIndex()] = flux;         // Put in flux array in main
                    if(showFluxCalc == 1){
                        s = "\n%d %18s reactants=%d iso0=%d iso1=%d iso2=%d Rrate=%7.3e Y1=%7.3e ";
                        s += "Y2=%7.3e Y3=%7.3e Flux=%7.3e";
                        printf(Utilities::stringToChar(s),
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], reactantIndex[1], 
                            reactantIndex[2], Rrate, Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], 
                            Y[ reactantIndex[2] ], flux);
                    }
                    break;
            }
            
        }    // End of function computeFlux()
        
        
        // Function Reaction::sumFplusFminus() to sum the total F+ and F- for each isotope.  
        
        static void sumFplusFminus(){
            
            //printf("\nCOMPUTING SUM F+ AND SUM F-\n");
            int minny = 0;
            double accum;
            double dydt;
            
            for(int i=0; i < numberSpecies; i++){	
                
                // Sum F+ for each isotope
                if(i>0) minny = FplusMax[i-1]+1;
                if(showFluxCalc == 1) printf("\n\nFplusMax=%d", FplusMax[i-1]);
                accum = zerod;	
                for(int j=minny; j<=FplusMax[i]; j++){
                    accum += Fplus[j];
                    if(showFluxCalc == 1) printf("\ni=%d j=%d Fplus=%g FplusSum=%8.4e", 
                        i, j, Fplus[j], accum);
                }
                setSpeciesfplus(i, accum);        // Also sets FplusSum[i] = accum;
                //printf("\nFplusSum[%d]=%8.4e", i, FplusSum[i]);
                
                // Sum F- for each isotope
                minny = 0;
                if(i>0) minny = FminusMax[i-1]+1;
                if(showFluxCalc == 1) printf("\n\nFminusMax=%d", FminusMax[i-1]);
                accum = zerod;
                for(int j=minny; j<=FminusMax[i]; j++){
                    accum += Fminus[j];
                    if(showFluxCalc == 1) printf("\ni=%d j=%d Fminus=%g FminusSum=%8.4e", 
                        i, j, Fminus[j], accum);
                }
                setSpeciesfminus(i, accum);      // Also sets FminusSum[i] = accum;
                setSpeciesdYdt(i, FplusSum[i] - FminusSum[i]);
                //printf("\nFminusSum[%d]=%8.4e", i, FminusSum[i]);
                
            }
            
            //printf("\n");
        }
        
        
};  // End class Reaction



// Class ReactionVector with  methods that create reaction vectors. Inherits from 
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
            printf("\n\nREACTION VECTOR ARRAY (%d Reaction vectors with %d species components):\n", 
                numberReactions, ISOTOPES);
            printf("\nReaction \\ Species           ");
            for(int k=0; k<uppity; k++){
                printf("%4s  ",isoLabel[k]);
            }
            printf("\n");
            for(int j=0; j<numberReactions; j++){
                printf("%4d %22s [ ",j,reacLabel[j]);
                for(int k=0; k<uppity-1; k++){
                    printf("%2d    ", reacMask[k][j]);
                }
                printf("%2d ]\n", reacMask[uppity-1][j]);
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
                printf("\nPopulate vector components of array rv[i]_j with reacMask[j][i]:");
            }
            
            for (int i = 0; i < SIZE; i++) {
                if(showRVdetails == 1) printf("\n\nrv[%d]",i);
                for(int j=0; j<ISOTOPES; j++){
                    gsl_vector_set (rvPt+i, j, reacMask[j][i]);
                    // Retrieve the vector component just stored and print it
                    int temp = gsl_vector_get(rvPt+i, j);
                    if(showRVdetails == 1){
                        printf("\ni=%d j=%d  reacMask[%d][%d] =%3d  rv[%d]_%d =%3d", 
                            i, j, i, j, reacMask[j][i], i, j, temp);
                    }
                }
            }
        
            
            // Display reaction vectors as component list
            
            printf("\nGSL REACTION VECTOR COMPONENTS (%d reaction vectors with %d components)\n",
                SIZE, ISOTOPES);
            
            for (int i = 0; i < SIZE; i++) {
                printf("\nrv[%d]: [",i);
                for(int j=0; j<ISOTOPES; j++){
                    
                    // Define a pointer that will point to the GSL vector in array entry rv[i].
                    gsl_vector *vector = rvPt+i;
                    
                    // Assign the jth component of the vector in rv[i] to a variable
                    int component = gsl_vector_get (vector, j);
                    printf ("%3d", component);
                }
                printf(" ]");
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
            
            if(showRGsorting == 1) printf("\n\n\n--- SORTING REACTION GROUPS ---");
            
            int scorekeeper = 0;
            for (int i=0; i<SIZE; i++){
                scorekeeper = 0;
                if(i==0) rindex ++;
                if(showRGsorting == 1) printf("\n\nRG=%d", rindex);
                for(int j=0; j<SIZE; j++){
                    
                    if(RGindex[i] < 0) RGindex[i] = rindex;
                    ck = compareGSLvectors(rvPt+i, rvPt+j);
                    
                    if(ck > 0 && RGindex[j]< 0) {
                        RGindex[j] = rindex;
                        scorekeeper ++;
                    }
                    if(showRGsorting==1){
                        printf("\ni=%d j=%d RGindex[%d]=%d ck=%d rindex=%d scorekeeper=%d", 
                            i, j, j, RGindex[j], ck, rindex, scorekeeper);
                    }
                }
                
                if(scorekeeper > 0) rindex++;
            }
            
            numberRG = rindex;   // Store total number of reaction groups
            
            // Diagnostic showing reaction group associated with each reaction
            
            if(showRGsorting == 1){
                printf("\n\n-- SUMMARY OF REACTION GROUPS:\n");
                for(int i=0; i<SIZE; i++){
                    printf("\nreaction=%d  %18s RGindex=%d RGmemberIndex=%d", 
                        i, reacLabel[i], RGindex[i], RGMemberIndex[i]);
                }
            }
            
            // Write out the components of the reaction groups
            
            printf("\n\n\nPARTIAL EQUILIBRIUM REACTION GROUPS");
            for(int i=0; i<numberRG; i++){
                printf("\n\nReaction Group %d:", i);
                int rgindex = -1;
                for(int j=0; j<SIZE; j++){
                    if(RGindex[j] == i){
                        rgindex ++; 
                        setRG(j, RGclass[j], RGindex[j]);
                        printf("\n%s reacIndex=%d RGindex=%d RG=%d RGreacIndex=%d isForward=%d RG: %s", 
                            reacLabel[j], j, rgindex, RGclass[j], RGMemberIndex[j],
                            isPEforward[j],
                            stringToChar(RGstring[j]));
                    }
                }
            }
            printf("\n");
            
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
                printf("\n\n--- Use parseF() to find F+ and F- flux components for each species ---");
            
            int incrementPlus = 0;
            int incrementMinus = 0;
            
            // Loop over all isotopes in the network		
            for(int i=0; i<numberSpecies; i++)
            {
                int total = 0;
                int numFplus = 0;
                int numFminus = 0;
                if(showParsing == 1) printf("\n");
                
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
                            printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F-", 
                                isoLabel[i], j, reacLabel[j], NumReactingSpecies[j], NumProducts[j], totalL, 
                                totalR, total);
                    } 
                    else if(total < 0){          // Contributes to F+ for this isotope
                        numFplus ++;
                        reacMask[i][j] = -total;
                        tempInt1[incrementPlus + numFplus-1] = j;
                        if(showParsing == 1)
                            printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F+", 
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
                    printf("\nSpecies=%d %s numF+ = %d numF- = %d", i, isoLabel[i], numFplus, numFminus);
            }
            
            // Display isotope component array
            
            printf("\n\n\nFLUX-ISOTOPE COMPONENT ARRAY (negative n for F-; positive n for F+ for given isotope):");
            printf("\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
            
            int uppity = minimumOf(30, numberSpecies);  // limit printout width to 30 species
            printf("\n\nIndex             Reaction");
            for(int k=0; k<uppity; k++){
                printf("%5s", isoLabel[k]);
            }
            for(int j=0; j<numberReactions; j++){
                printf("\n%3d %22s",j,reacLabel[j]);
                for(int k=0; k<uppity; k++){
                    printf(" %4d",reacMask[k][j]);
                }
            }
            
            printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.\n", 
                totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
            
        }   // End of function parseF()
    
    
};    // end class ReactionVector


/*
*   class MatrixUtils inherits from Utilities. Methods to create GSL matrices and vectors
*   given standard C++ arrays and to compute matrix-vector multiply M*v using BLAS.
*/


class MatrixUtils: public Utilities {

    private:

        // Rows and columns in flux matrix
        // Number of columns in flux matrix is always to rows in abundance vector
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
        int niso;                              // Number of isotopic species in RG class, this object
        int RGn;                               // Index of reaction group in RG array (0, 1, ... #RG)
        int numberMemberReactions;             // Number of reactions in this RG instance
        int memberReactions[maxreac];          // reacIndex of reactions in reaction group
        int numberReactants[maxreac];          // Number of reactants for each reaction in RG
        int numberProducts[maxreac];           // Number of products for each reaction in RG

        int rgclass;                           // Reaction group class (0-5)
        bool isEquil;                          // True if RG in equilibrium; false otherwise
        bool isEquilMaybe;                     // Whether would be in equil if no threshhold condition
        bool isForward[maxreac];               // Whether reaction in RG labeled forward
        double flux[maxreac];                  // Current flux for each reaction in RG
        double netflux;                        // Net flux for the entire reaction group
        char reaclabel[maxreac][LABELSIZE];    // Member reaction label
        int RGarrayIndex;                      // Index of ReactionGroup RG[] array
        
        // Partial equilibrium quantities
        
        double Yzero[ISOTOPES];        // Hold Y for this species at beginning of timestep
        
        double crg[4];                 // Constants c1, c2, ... (1-4 entries; could allocate dynamically)
        int numberC;                   // Number of constants crg[] for this rg class (1-4 entries)
        double rgkf;                   // Forward rate parameter for partial equilibrium
        double rgkr;                   // Reverse rate parameter for partial equilibrium
        
        double aa, bb, cc;             // Quadratic coefficients a, b, c
        double alpha, beta, gamma;     // Helper coefficients for cubic ~ quadratic approximation
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
        int isoZ[5];                   // Z for the niso isotopes in the reactions of the group
        int isoN[5];                   // N for the niso isotopes in the reactions of the group
        double isoA[5];                // A for the niso isotopes in the reactions of the group
        double isoYeq[5];              // Y_eq for the niso isotopes in the reactions of the group
        double isoY[5];                // Current Y for the niso isotopes in the reactions of the group
        double isoY0[5];               // Current Y for isotopes in the reactions of the group

    
    public:
    
    // Constructor
    
    ReactionGroup(int rgn){
        RGn = rgn;
        isEquil = false;
        isEquilMaybe = false;
    }
    
    // Public ReactionGroup setter functions to set values of private class fields
    
    void setnumberMemberReactions(int n){numberMemberReactions = n;}
    
    void setmemberReactions (int i, int index){memberReactions[i] = index;}
    
    void setnumberReactants(int i, int j){numberReactants[i] = j;}
    
    void setnumberProducts(int i, int j){numberProducts[i] = j;}
    
    void setisEquil(bool b) {isEquil = b;}
    
    void setisForward(int i, bool b){isForward[i] = b;}
    
    void setflux(int i, double f){flux[i] = f;}
    
    void setRGn(int rgn){RGn = rgn;}
    
    void setrgclass(int rc){rgclass = rc;}

    void setreacString(int k, string s){ 
        // Convert from string to char array
        char p[s.length()+1];  
        for (int i = 0; i < sizeof(p); i++) { 
            p[i] = s[i]; 
            reaclabel[k][i] = p[i];
        }
    }
       
    // Method to set all fluxes in RG
    
    void setRGfluxes(){
        //printf("\n\n**** setRGFluxes() t = %7.4e", t);
        for(int i=0; i<numberMemberReactions; i++){
            setflux(i, Flux[ memberReactions[i] ]);
        }
    }
    
    void setnetflux(double f){netflux = f;}
    
    void setRGarrayIndex(int i){RGarrayIndex = i;}
    
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
    
    double getflux(int i){return flux[i];}
    
    int getRGn(){return RGn;}
    
    int getrgclass(){return rgclass;}
    
    char* getreacString(int k){return reaclabel[k];}
    
    double getnetflux(){return netflux;}
    
    int getRGarrayIndex(){return RGarrayIndex;}
    
    int getniso(){return niso;};
    
    int getisoindex(int i){return isoindex[i];}
    
    char* getisolabel(int k){return isolabel[k];}
    
    int getisoZ(int i){return isoZ[i];}
    
    int getisoN(int i){return isoN[i];}
    
    int getisoA(int i){return isoA[i];}
    
    int getreactantIsoIndex(int i){return reactantIsoIndex[i];}
    
    int getproductIsoIndex(int i){return productIsoIndex[i];}
    
    double geteqcheck(int k){return eqcheck[k];}
    
    
    // Method to show all current fluxes in reaction groups
    
    void showRGfluxes(){
        
        //printf("\n**** showRGFluxes() t = %7.4e\n", t);
        
        double fac;
        for(int i=0; i<numberMemberReactions; i++){
            if(getisForward(i)){
                fac = 1.0;
            } else {
                fac = -1.0;
            }
//             printf("*****memberIndex=%d %s RGclass=%d isForward=%d flux=%7.4e\n", 
//                 i, 
//                 reacLabel[memberReactions[i]],
//                 getrgclass(),
//                 getisForward(i),   // prints 1 if true; 0 if false
//                 fac*getflux(i)
//             );
        }
        if(isEquil){
            printf("RG=%d NetRGflux=%7.4e (Equilibrated)\n", RGn, netflux); 
        } else {
            printf("RG=%d NetRGflux=%7.4e (Not Equilibrated)\n", RGn, netflux); 
        }
    }
    
    // Method to sum net flux for this reaction group
    
    double sumRGfluxes(){
        
        //printf("\n**** sumRGFluxes() t = %7.4e\n", t);
        
        double sumf = zerod;
        double fac;

        for (int i=0; i<numberMemberReactions; i++){
            fac = -1.0;
            if(isForward[i]) fac = 1.0;
            sumf += fac*flux[i];
        }
        //printf("\n");
        netflux = sumf;
        return sumf;
        
    }    // End function sumRGfluxes()
    
    
    // ---------------------------------------------------------
    // Method to compute all partial equilibrium quantities
    //----------------------------------------------------------
    
    void computeEquilibrium() {
        
        //printf("\n****Compute Equilibrium");
        
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
    // Method to compute the net forward and reverse rates k_f and k_r
    // required in partial equilibrium approximation.
    // -----------------------------------------------------------------
    
    void computeEquilibriumRates() {
        
        //printf("\n****computeEquilibriumRates()");
        
        double kf = 0;
        double kr = 0;
        
        // Sum all contributions from members of reaction group
        
        for (int j = 0; j < numberMemberReactions; j++) {
            int rin = memberReactions[j];
            if (isForward[j]) {
                kf += Rate[rin];
//                 kf += masterRates[this.reactions[j].Z][this.reactions[j].N]
//                 [this.reactions[j].reacIndex];
            } else {
                kr += Rate[rin];
//                 kr += masterRates[this.reactions[j].Z][this.reactions[j].N]
//                 [this.reactions[j].reacIndex];
            }
        }
        
        // Store forward and reverse rates
        
        rgkf = kf;
        rgkr = kr;
        
    }  // End of function computeEquilibriumRates()
    
    
    // -----------------------------------------------------------------------
    // Method to put the values of Y0 at beginning of timestep into the Y0[]
    // array for this object
    // -----------------------------------------------------------------------
    
    void putY0() {
        
        printf("\n**** Put Y0 niso=%d RGarrayindex=%d\n", niso, RGarrayIndex);
        
        int ii;
        
        for (int k = 0; k < niso; k++) {    // loop over the niso isotopes in RG
            ii = isoindex[k];
            isoY0[k] = Y[ii];
            isoY[k] = isoY0[k];
            printf("\nRG=%d k=%d isoindex=%d isoY0[%s]=%7.3e", 
                   RGarrayIndex, k, ii, isoLabel[ii],  isoY[k]);
        }
        printf("\n");
    }
    
   
   
    // Method to compute the values of the constants crg[]
    
    void computeC() {
        
        printf("\n****computeC()");
        
        switch (rgclass) {
            
            // Reaclib class 7, which can't equilibrate
            case -1: 
                break;
                
            case 1:
                crg[0] = isoY0[0] + isoY0[1];
                break;
                
            case 2:
                crg[0] = isoY0[1] - isoY0[0];
                crg[1] = isoY0[1] + isoY0[2];
                break;
                
            case 3:
                crg[0] = isoY0[0] - isoY0[1];
                crg[1] = isoY0[0] - isoY0[2];
                crg[2] = THIRD * (isoY0[0] + isoY0[1] + isoY0[2]) + isoY0[3];
                break;
                
            case 4:
                crg[0] = isoY0[0] - isoY0[1];
                crg[1] = isoY0[0] + isoY0[2];
                crg[2] = isoY0[0] + isoY0[3];
                break;
                
            case 5:
                crg[0] = isoY0[0] + THIRD * (isoY0[2] + isoY0[3] + isoY0[4]);
                crg[1] = isoY0[0] - isoY0[1];
                crg[2] = isoY0[2] - isoY0[3];
                crg[3] = isoY0[2] - isoY0[4];
                break;
        }
    }     // End of method computeC()
    
    
    // --------------------------------------------------------------
    // Method ReactionGroup::computeQuad() to compute the quadratic 
    // coefficients needed for the equilibrium solution and to compute 
    // the equilibrium solution
    // --------------------------------------------------------------
    
    void computeQuad() {
        
        printf("\n****computeQuad()");
        
        switch (rgclass) {
            
            // Reaclib class 7, which can't equilibrate
            case -1: 
                
                break;
                
            case 1:
                aa = 0;
                bb = -rgkf;
                cc = rgkr;
                break;
                
            case 2:
                aa = -rgkf;
                bb = -(crg[0] * rgkf + rgkr);
                cc = rgkr * (crg[1] - crg[0]);
                break;
                
            case 3:
                aa = -rgkf * isoY0[0] + rgkf * (crg[0] + crg[1]);
                bb = -(rgkf * crg[0] * crg[1] + rgkr);
                cc = rgkr * (crg[2] + THIRD * (crg[0] + crg[1]));
                break;
                
            case 4:
                aa = rgkr - rgkf;
                bb = -rgkr * (crg[1] + crg[2]) + rgkf * crg[0];
                cc = rgkr * crg[1] * crg[2];
                
                break;
                
            case 5:
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
        } else {
            qq = -1;
            tau = 1 / rgkf;
            isoYeq[0] = rgkr / rgkf;
        }
        
        // Compute the other equilibrium populations in the reaction pair
        // and abundance ratios
        
        switch (rgclass) {
            
            // Reaclib class 7, which can't equilibrate
            case -1: 
                
                break;
                
            case 1:
                isoYeq[1] = crg[0] - isoYeq[0];
                equilRatio = isoY[0] / isoY[1];
                break;
                
            case 2:
                isoYeq[1] = crg[0] + isoYeq[0];
                isoYeq[2] = crg[1] - isoYeq[1];
                equilRatio = isoY[0] * isoY[1] / isoY[2];
                break;
                
            case 3:
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = isoYeq[0] - crg[1];
                isoYeq[3] = crg[2] - isoYeq[0] + THIRD * (crg[0] + crg[1]);
                equilRatio = isoY[0] * isoY[1] * isoY[2] / isoY[3];
                break;
                
            case 4:
                isoYeq[1] = isoYeq[0] - crg[0];
                isoYeq[2] = crg[1] - isoYeq[0];
                isoYeq[3] = crg[2] - isoYeq[0];
                equilRatio = isoY[0] * isoY[1] / (isoY[2] * isoY[3]);
                break;
                
            case 5:
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
        
    }    // End method computeQuad()
    
    
    // ---------------------------------------------------------------------
    // Method to compute q = 4ac-b^2 for quadratic solution
    // ---------------------------------------------------------------------
    
    double computeq(double a, double b, double c) {
        return 4 * a * c - b * b;
    }
    
    // Method to compute Yeq[0]
    
    double computeYeq(double a, double b, double rootq) {
        return -0.5 * (b + rootq) / a;
    }
    
    
    // ---------------------------------------------------------------------
    // Method to compute array of population ratios used to check
    // equilibration
    // ---------------------------------------------------------------------
    
    void computeEqRatios() {
        
        printf("\n****computeEqRatios()\n");
        
        double thisDevious = abs((equilRatio - kratio) / kratio);
        
        if (isEquil && thisDevious > mostDevious) {
            mostDevious = thisDevious;
            mostDeviousIndex = RGarrayIndex;
        }
        
        // The return statements in the following if-clauses cause reaction
        // groups already in equilibrium to stay in equilibrium. If the 
        // maxDevious > tolerance check is implemented it can cause a
        // reaction group to drop out of equilibrium.
        
        if (isEquil && thisDevious < maxDevious) {
            return;
        } else if (isEquil && thisDevious > maxDevious && imposeEquil
            && t > equilibrateTime) {
            removeFromEquilibrium();
            return;
        }
            
            Yminner = 1000;
            maxeqcheck = 0;
            mineqcheck = 1000;
            
            // Determine if RG is in equilibrium: set isEquil to default value
            // of true and then try to falsify
            
            isEquil = true;
            isEquilMaybe = true;
            
            for (int i = 0; i < niso; i++) {
                
                // Note: something like the following probably required because
                // otherwise we will divide by zero for isotopes early in the 
                // calculation that have no population.
                
                if (isoYeq[i] == 0 || isoY[i] == 0) {
                    isEquil = false;
                    isEquilMaybe = false;
                    break;
                }
                
                eqcheck[i] = abs(isoY[i] - isoYeq[i]) / isoYeq[i];
                
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
                
                if (t < equilibrateTime || eqcheck[i] > equiTol
                    || isoYeq[i] < Ythresh) {
                    isEquil = false;
                
                    // break; // Note: this break won't affect results, but
                    // would affect diagnostic values of eqcheck[]
                    // since they will all be zero after the break.
                
                }
            }
            
            // Check whether would be in equil without time or threshhold condition
            
            for (int i = 0; i < niso; i++) {
                
                if (eqcheck[i] > equiTol) {
                    isEquilMaybe = false;
                    break;
                }
            }
            
            /* Keep track of number of reaction in partial equilibrium.  totalEquilReactions
             *  is only updated here if imposeEquil is false. If imposeEquil is true,
             *  totalEquilReaction is updated when the flux is suppressed for equilibrium 
             *  pairs. */
            
            if (!imposeEquil && isEquil)
                totalEquilReactions += numberMemberReactions;
                totalEquilRG ++;
            
            // Set isEquil field of network species vectors to true if isotope
            // participating in equilibrium
            
            if (isEquil) {
                if (showAddRemove) {
                    printf("\n*** totalTimeSteps=%d Adding RG %d to equilibrium\n",
                        totalTimeSteps, RGarrayIndex);
                }
                for (int i = 0; i < niso; i++) {
                    //netVector[this.abundVecIndex[i]].isEquil = true;
                }
            }
            
            // Set the activity array for each reaction in reaction group to true if not in 
            // equil and false if it is, if we are imposing equilibrium.
            
            if (imposeEquil && t > equilibrateTime) {
                for (int i = 0; i < numberMemberReactions; i++) {
                    int ck = memberReactions[i];
                    reacIsActive[ck] = !isEquil;
                }
            }
            
    }    // End method computeEqRatios()
    
    
    // -----------------------------------------------------------
    // Method to remove reaction group from equilibrium
    // -----------------------------------------------------------
    
    void removeFromEquilibrium() {
        
        printf("\n ****removeFromEquilibrium()");
        
        isEquil = false;
        double thisDevious = abs((equilRatio - kratio) / kratio);
        if (showAddRemove) {
            printf("\n*** %d Remove RG %d devious=%d Rmin=%8.4e Rmax=%8.4e Ymin=%8.4e\n", 
                RGarrayIndex, thisDevious, mineqcheck, maxeqcheck, Yminner);
        }
        for (int i = 0; i < niso; i++) {
            isEquil = false;
            if (showAddRemove) {
                printf("\nZ=%d N=%d Y=%8.4e Yeq=%8.4e Rprev=%8.4e Rnow=%8.5e",
                    isoZ[i], isoN[i], isoY[i], isoYeq[i], eqcheck[i],
                    abs(isoY[i] - isoYeq[i]) / isoYeq[i]
                );
            }
        }
        
        for (int i = 0; i < numberMemberReactions; i++) {
            int ck = memberReactions[i];
            reacIsActive[ck] = true;         
            if (showAddRemove) {
                printf("\n *** Remove %s", reacLabel[i]);
            }
        }
    }
    
    
    // ----------------------------------------------------------------
    // Method to determine if given species with index speciesIndex is 
    // in any of the reactions of a reaction group, where speciesIndex
    // is the array index i for isotope quantitites like Z[i]. Returns 
    // true (1) if it is and false (0) if not.
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



// Class Integrate with methods to integrate the reaction network.  Inherits from
// class Utilities.

class Integrate: public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions.  Its static functions can be called directly from the class
    // without having to instantiate.
    
    private:
        
    public:
        
        // Function to execute a single integration step.  Assumes that all fluxes have
        // already been calculated and relevant fluxes set to zero if PE approximation
        // and the reaction group has been judged to be in equilibrium.
        
        static void doIntegrationStep(){
            
            // Determine trial timestep.  This timestep will be updated in the various 
            // integration methods to a final timestep using getTimestep(). 
            
            if(constantTimestep){
                dt = constant_dt;      // Constant timestep
            } else {
                dt = 0.01 * t;     // Temporary for testing
                //dt = getTrialTimestep();    // Trial adaptive timestep
                //dt = getTimestep();    // Adaptive timestep
            }
            
            // If using the QSS approximation, apply QSS approximation to all isotopes
            
            if(doQSS){
                printf("\nQSS approximation (t=%7.4e, dt=%7.4e)\n", t, dt);
                for (int i=0; i<ISOTOPES; i++){
                    QSSupdate(Fplus[i], Fminus[i], Y[i], dt);
                    printf("\n-------Doing QSS update for %s  Y=%7.4e", isoLabel[i], Y[i]);
                }
                printf("\n");
                
            }
            
            // If using asymptotic approximation, determine which species satisfy the
            // asymptotic condition.
            
            if(doASY){
        
//                 printf("\nCheck asymptotic condition (t=%7.4e, dt=%7.4e)\n",
//                     t, dt
//                 );
                for(int i=0; i<ISOTOPES; i++){
                    isAsy[i] = checkAsy(FminusSum[i], Y[i], dt);
                    if(isAsy[i]){ totalAsy++; } 
                }
                
                // Summarize results
                
                if(showAsyTest){
                    //printf("\nindex   iso   FminusSum           Y       check    Asy");
                    string asyck;
                    for(int i=0; i<ISOTOPES; i++){
                        asyck = "false";
                        if(isAsy[i]) asyck="true";
                        double ck;
                        if(Y[i] > zerod){
                            ck = FminusSum[i]*dt/Y[i];
                        } else {
                            ck = zerod;
                        }
                        
//                         printf("\n    %d %5s  %7.4e  %7.4e  %7.4e  %5s",
//                             i, (isoLabel[i]), FminusSum[i], 
//                                Y[i], ck, Utilities::stringToChar(asyck)
//                         );
                    }
                    //printf("\n");
                }
                
                // If Asy+PE, compute the matrix multiply for forward euler, with fluxes removed
                // by PE approximation (individual fluxes in RG that are in equilibrium) and 
                // asymptotic approximation (rows and columns of matrix)
                
                updateAsyEuler();
            }
            
            
            
            
            if(!doQSS){
                // Call matrix multiply with elements of matrix having been removed by
                // PE and Asy approximations.
            }
            
            // Check the sum of the mass fractions. Should be 1.0 if particle number is
            // being conserved
            
            sumX = Utilities::returnSumX();
            diffX = sumX - unitd;
            
        }    // end of doIntegrate(I)
        
        

        
        // Function to set the trial integration timestep
        
        static double getTrialTimestep(){
            
            // Placeholder. This will supply the initial trial timestep,
            // which is required to initiate the iteration to the final 
            // timestep for this time interval.
            
            double timestep = constant_dt;
            return timestep;
        }
        
        
    // Function to set the current integration timestep
        
    static double getTimestep(){
        
        // Placeholder. The adaptive timestepping algorithm will go here. We should be
        // able to use the neutrino transport adaptive timestepper that Aaron and Adam
        // have been testing with suitable modification.  For now just 
        // return a constant.
        
        double timestep = constant_dt;
        return timestep;
    }
    
     
    // Function to update by the forward Euler method
        
    static double eulerUpdate(double FplusSum, double FminusSum, double Y, double dt){
        
//         printf("Forward Euler input: FplusSum = %9.5e FminusSum = %9.5e Y = %9.5e dt = %9.5e\n", 
//                FplusSum, FminusSum, Y, dt);
        return Y + (FplusSum-FminusSum)*dt;   // New Y for forward Euler method
        
    }
    
    // Function to update by the asymptotic method
    
    static double asymptoticUpdate(double Fplus, double Fminus, double Y, double dt){
        
        // Update Y by asymptotic approximation (Sophia He formula)
        
        printf("Asymptotic input: Fplus = %9.5e Fminus = %9.5e Y = %9.5e dt = %9.5e", 
            Fplus, Fminus, Y, dt);
        
        return (Y + Fplus*dt)/(unitd + Fminus*dt/Y);  // New Y for asymptotic method
        
    }
    
    // Function to determine whether an isotope satisfies the
    // asymptotic condition. Returns true (1) if it does and false (0) if not.
    
    static bool checkAsy(double Fminus, double Y, double dt){
        
//         printf("Asymptotic check input: Fminus = %9.5e Y = %9.5e dt = %9.5e ck=%9.5e", 
//             Fminus, Y, dt, Fminus*dt/Y);
        
        if(Y > zerod && Fminus*dt/Y > unitd){
            return true;
        } else {
            return false;
        }
        
    }
    
    
    // Function to update by the Quasi-Steady-State (QSS) approximation.  Placeholder
    // for now.  Should be able to adapt Adam's neutrino QSS algorithm.
    
    static void QSSupdate(double Fplus, double Fminus, double Y, double dt){
        
        // *************************
        // QSS algorithm goes here.
        // *************************
        
    }
    
    
    /* Use the fluxes to update the populations for this timestep
     For now we shall assume the asymptotic method. We determine whether each isotope 
     satisfies the asymptotic condition. If it does we update with the asymptotic formula. 
     If not, we update numerically using the forward Euler formula. */
    
    static void updateAsyEuler(){
        //printf("\n\n$$$$$ Updating asy-euler\n");
        for(int i=0; i<numberSpecies; i++){		
            if(checkAsy(Fminus[i], Y[i], dt) == 1){
                Y[i] = asymptoticUpdate(FplusSum[i], FminusSum[i], Y[i], dt);
            } else {
                Y[i] = eulerUpdate(FplusSum[i], FminusSum[i], Y[i], dt);
            }
            X[i] = Y[i] * (double) AA[i];
        }
    }    // End function updateAsyEuler()
    
    
};    // End class integrate

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
    
    // Set labels and check consistency of choice for explicit algebraic methods set.
    // Generally we use either asymptotic (Asy) or quasi-steady-state (QSS) algorithms.
    // In either case we may choose to add the partial equilibrium (PE) algorithm. So
    // valid options are Asy, QSS, Asy+PE, and QSS+PE.
    
    methstring="Using ";
    if(doASY){
       methstring += "ASY";
       doQSS = false;
    } else {
        methstring += "QSS";
        doASY = false;
    }
    if(doPE){methstring += "+PE";}
    methstring += " method";
    printf("%s\n", Utilities::stringToChar(methstring));
    
    // Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
    // realistic calculation the temperature and density will be passed from the hydro 
    // code in an operator-split coupling of this network to hydro. Here we hardwire
    // them for testing purposes.  These will be used to calculate the reaction
    // rates in the network. If we assume operator splitting, the temperature 
    // and density are assumed constant for each network integration. But we allow
    // the possibility below to interpolate the temperature and density from a
    // hydrodynamical profile as a function of time.
    
    T9_start = 5.0f;
    T9 = T9_start;
    rho_start = 1.0e8;
    rho = rho_start;
    
    // Set the range of time integration and the initial timestep (units of seconds).  
    // In an operator-split coupling tmax will come from the hydro and dt_init will 
    // likely be the last timestep of the previous network integration (for the preceding 
    // hydro timestep). Here we hardwire them for testing purposes.
    
    double tmax = 1e-11;
    double dt_init = 1e-17; 
    
    // Read in network file and associated partition functions.  This is required only
    // once at the beginning of the entire calculation.  
    
    char *networkFilePtr = networkFile;
    readNetwork(networkFilePtr);
    writeNetwork();
    
    // Read in rate library data from a file. This is required only once, at the
    // beginning of the entire calculation.
    
    char *rateLibraryFilePtr = rateLibraryFile;
    readLibraryParams(rateLibraryFilePtr);
    
//     // Print out some quantitites from the Reaction object reaction[].  
    
    for(int i=0; i<SIZE; i++){
        
        printf("\n%d %s reacClass=%d reactants=%d products=%d isEC=%d isReverse=%d Q=%5.4f prefac=%5.4f", 
            reaction[i].getreacIndex(), 
            Utilities::stringToChar(reaction[i].getreacString()),  
            reaction[i].getreacClass(),
            reaction[i].getnumberReactants(),
            reaction[i].getnumberProducts(),
            reaction[i].getisEC(),
            reaction[i].getisReverse(),
            reaction[i].getQ(),
            reaction[i].getprefac()
        );
    }
    
    printf("\n\n\nREACLIB PARAMETERS FOR %d REACTIONS\n", SIZE);
    printf("\n                                p0         p1         p2         p3         ");
    printf("p4         p5         p6");
    for (int i=0; i<SIZE; i++){
        printf("\n%3d  %18s %10.4f", i, reaction[i].getreacChar(), reaction[i].getp(0));
        for(int j=1; j<7; j++){
            printf(" %10.4f", reaction[i].getp(j));
        }
    }
    
    printf("\n\nZ and N for reactants:\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        printf("\n%3d %18s ", i, reaction[i].getreacChar());
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            printf(" Z[%d]=%d", j, reaction[i].getreactantZ(j));
        }
        printf(" ");
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            printf(" N[%d]=%d", j, reaction[i].getreactantN(j));
        }
    }
    
    printf("\n\nZ and N for products:\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        printf("\n%3d %18s ", i, reaction[i].getreacChar());
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            printf(" Z[%d]=%d", j, reaction[i].getproductZ(j));
        }
        printf(" ");
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            printf(" N[%d]=%d", j, reaction[i].getproductN(j));
        }
    }

    printf("\n\nreactantIndex for %d reactions (index of species vector for each reactant):\n", SIZE);

    for (int i=0; i<SIZE; i++){
        printf("\n%d %18s ",i,reaction[i].getreacChar());
        for(int j=0; j<reaction[i].getnumberReactants(); j++){
            printf(" reactantIndex[%d]=%d", j, reaction[i].getreactantIndex(j));
        }
    }
    
    printf("\n\nproductIndex for %d reactions (index of species vector for each product):\n", SIZE);
    
    for (int i=0; i<SIZE; i++){
        printf("\n%d %18s ",i,reaction[i].getreacChar());
        for(int j=0; j<reaction[i].getnumberProducts(); j++){
            printf(" productIndex[%d]=%d", j, reaction[i].getproductIndex(j));
        }
    }
    
    // Find the time intervals for plot output during the integration. After this
    // function is executed the plotSteps target time intervals for output will
    // be in the array plotTimeTargets[]. In the integration the ith output step will 
    // be triggered as soon as the time t is >= plottimeTargets[i].  The actual time of the output
    // (which will usually be slightly larger than plottimeTargets[i]) will be stored in tplot[i].
    
    Utilities::log10Spacing(start_time, stop_time, plotSteps, plotTimeTargets);
    
    printf("\n\nPlot Intervals:\n");
    for(int i=0; i <plotSteps; i++){
        printf("\ni=%d tplot=%7.4e log(tplot)=%7.4e", 
               i, plotTimeTargets[i], log10(plotTimeTargets[i]));
    }
    
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
    printf("\n\nNETWORK SPECIES VECTOR (%d components):\n\nIndex  Species    Z     N",
        numberSpecies);
    for(int i=0; i<numberSpecies; i++){
        printf("\n%5d    %5s  %3d  %4d", i, isoLabel[i], Z[i], N[i]);
    }
    printf("\n");
    
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
    
    // Assign values of fields for the ReactionGroup objects RG[]
    
    assignRG();
    
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
    // Method declared static so it can be called as Reaction::setupFplusFminus() 
    // without having to instantiate.
    
    Reaction::setupFplusFminus();
    
    
    
    // -----------------------------------------------
    // *** Begin main time integration while-loop ***
    // -----------------------------------------------
    
    
    printf("\n\n\n                 --- BEGIN TIME INTEGRATION ---\n");
    
    dt = dt_start;              // Integration start time
    t = start_time - dt;        // Current integration time
    totalTimeSteps = 0;         // Integration step counter
    totalEquilRG = 0;           // Number quilibrated reaction groups
    totalEquilReactions = 0;    // Number equilibrated reactions
    totalAsy = 0;               // Number asymptotic species
    int plotCounter = 1;        // Plot output counter
    
    Utilities::startTimer();    // Start a timer for integration
    
    // Compute initial rates. If constant_T9 and constant_rho are true, rates won't
    // change in the integration and don't need to be computed again.  If either
    // T9 or rho change, the rates will be recomputed at each integration step.
    // Use methods of Reaction class to compute reaction rates. We have instantiated
    // a set of Reaction objects in the array reaction[i], one entry for each
    // reaction in the network. Loop over this array and call the computeRate()
    // method of Reaction on each object. 
    
    printf("\n\nINITIAL COMPUTED RATES\n");
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeConstantFacs(T9, rho);
        reaction[i].computeRate(T9, rho);
    }
    
    if(constant_T9 && constant_rho){
        printf("\n\n**** Rates not computed again since T and rho won't change in integration ****\n");
    } else {
        printf("\n\n**** Rates will be recomputed at each timestep since T and rho may change ****\n");
    }
    
    
    
    while(t < stop_time){
        
        t += dt;                
        totalTimeSteps ++;  
        
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
    
        // Use methods of Reaction class to compute reaction rates. We have instantiated
        // a set of Reaction objects in the array reaction[i], one entry for each
        // reaction in the network. Loop over this array and call the computeRate()
        // method of Reaction on each object. If constant_T9 and constant_rho are true, 
        // the rates only need be computed once as they won't change. Otherwise they are
        // recomputed for each integration step.
        
        if( (!constant_T9 || !constant_rho) && totalTimeSteps > 1){
            
            printf("**** RECOMPUTED RATES, timestep=%d\n",totalTimeSteps);
            
            for(int i=0; i<SIZE; i++){
                reaction[i].computeConstantFacs(T9, rho);
                reaction[i].computeRate(T9, rho);
            }
        }
        
        
        // Use methods of the Reaction class to compute fluxes.  We have instantiated
        // a set of Reaction objects in the array reaction[i], one entry for each
        // reaction in the network. Loop over this array and call the computeFlux()
        // method of Reaction on each object. Fluxes must be recomputed at each timestep
        // since they depend on the rates and the abundances. If temperature and density
        // are constant the rates won't change, but the fluxes generally will since
        // the abundances change even if the rates are constant.
        
        //printf("\n\nTOTAL FLUXES\n");
        
        for(int i=0; i<SIZE; i++){
            reaction[i].computeFlux();
        }
        
        //printf("\n\n\nPOPULATE REACTION GROUPS WITH FLUXES AND ABUNDANCES\n");
        
        for(int i=0; i<numberRG; i++){
            //printf("\nRG=%d", i);
            RG[i].setRGfluxes();
            if(doPE){
                RG[i].sumRGfluxes();
                RG[i].showRGfluxes();
                RG[i].computeEquilibrium();
            }
        }
        
        // If partial equilibrium approximation (doPE = true), set fluxes identically to
        // zero for all reactions in reaction groups that are judged to be in equilibrium
        // (RG[i].isEquil = true).
        
        if(doPE){
            
            //printf("\n\nIMPOSE EQUILIBRIUM CONDITION ON FLUXES");
            
            // Loop over reaction groups
            for(int i=0; i<numberRG; i++){
                
                // If RG equilibrated, loop over members of reaction group and set
                // each reaction flux to zero.
                
                bool ckequil = RG[i].getisEquil();
                //printf("\n\nRG=%d isEquil=%d", i, ckequil);
                
                if(ckequil){
                    for(int j=0; j<RG[i].getnumberMemberReactions(); j++){
                        Flux[RG[i].getmemberReactions(j)] = zerod;  // Set identically zero
                    } 
                }
                
                // Print results 
                
                for(int j=0; j<RG[i].getnumberMemberReactions(); j++){
//                     printf("\n   %d reaction=%d %s flux=%7.4e eqcheck=%7.4e",
//                         j, RG[i].getmemberReactions(j), 
//                         RG[i].getreacString(j),
//                         Flux[RG[i].getmemberReactions(j)],
//                          RG[i].geteqcheck(j)
//                     );
                }
            }
        }
        
        
        // Call the static function Reaction::populateFplusFminus() to populate F+ and F-
        // for each isotope set up in setupFplusFminus() from master flux array computed
        // with setRGfluxes() above.
        
        Reaction::populateFplusFminus();
        
        // Sum F+ and F- for each isotope
        
        Reaction::sumFplusFminus();
        
        // Perform an integration step
        
        Integrate::doIntegrationStep();
        
        // Display and output to files updated quantities at plotSteps times corresponding
        // to (approximately) equally-spaced intervals in log_10(time). The target output 
        // times are generated by Utilities::log10Spacing() and stored in the array
        // plotTimeTargets[plotSteps]. A screen and plot file output is triggered if
        // t >= plotTimeTargets[plotCounter-1], so the actual output times may be slightly
        // larger than the target output times.
        
        if(t >= plotTimeTargets[plotCounter-1]){
            
            string dasher2 = dasher + dasher;   // Dashed-line separator
            
            // Output to screen
            
            printf("\n%s", Utilities::stringToChar(dasher2));
            printf("\n%d/%d steps=%d T9=%4.2f rho=%4.2e t=%8.4e log_t=%6.4f dt=%8.4e sumX=%6.4f", 
                plotCounter, plotSteps, totalTimeSteps, T9, rho, t, log10(t), dt, sumX);
            printf("\n%s", Utilities::stringToChar(dasher2));
            tempest = "\nIndex   Iso           Y           X        dY/dt";
            tempest += "        dX/dt           dY           dX\n";
            printf(Utilities::stringToChar(tempest));
            
            for(int i=0; i<ISOTOPES; i++){
                printf("%5d %5s  %8.4e  %8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e\n", 
                    i, isoLabel[i], Y[i], X[i], isotope[i].getdYdt(), isotope[i].getdXdt(),
                       isotope[i].getdYdt()*dt, isotope[i].getdXdt()*dt
                );
            }
            printf("%s\n", Utilities::stringToChar(dasher2));
            
            // Output to plot arrays for this timestep
            
            tplot[plotCounter-1] = log10(t);
            dtplot[plotCounter-1] = log10(dt);
            
            // Following 2 temporary placeholders until energy calculations inserted
            EReleasePlot[plotCounter-1] = 20.0; //log10(abs(ERelease));
            dEReleasePlot[plotCounter-1] = 20.0; //log10(abs(dERelease));
            
            sumXplot[plotCounter-1] = sumX;
            numAsyplot[plotCounter-1] = totalAsy;
            totalEquilRG = 0;
            for(int i=0; i<numberRG;i++){
                if(RG[i].getisEquil()) totalEquilRG ++;
            }
            numRG_PEplot[plotCounter-1] = totalEquilRG;
            
            for(int i=0; i<ISOTOPES; i++){
                Xplot[i][plotCounter-1] = X[i];
            }
            
            //printf("\nplotSteps=%d logt=%6.3f\n", plotCounter, tplot[plotCounter]);
            
            // Increment the plot counter for next output
            plotCounter ++;
        }
    
    }   // End time integration while-loop
    
    
    printf("\nEnd of integration");
    Utilities::stopTimer();        // Stop timer and print time for integration
    printf("\n");

    // ------------------------------
    // *** End time integration ***
    // ------------------------------
    

    // Display abundances and mass fractions at end of integration

    printf("\nFINAL ABUNDANCES Y AND MASS FRACTIONS X\n");

    for(int i=0; i<ISOTOPES; i++){
        tempest = "\n%d %s Y=%7.3e X=%7.3e F+Sum=%7.3e ";
        tempest += "F-Sum=%7.3e dY/dt=%+7.3e dX/dt=%+7.3e";
        printf(Utilities::stringToChar(tempest), 
               i, 
               isotope[i].getLabel(), 
               isotope[i].getY(), 
               isotope[i].getX(),
               isotope[i].getfplus(),     // or FplusSum[i],
               isotope[i].getfminus(),    // or FminusSum[i]
               isotope[i].getdYdt(),
               isotope[i].getdXdt()
        );
    }

    printf("\n\n");
    
    // Output of data to plot files after integration
    
    Utilities::plotOutput();
    
    
    // ------------------------------------------------------------
    // --- Perform some optional tests of various functions ---
    // ------------------------------------------------------------
    
    if(showFunctionTests==1){
    
        printf("\n\n-- TEST OF SOME FUNCTIONS --\n");
        
        // Demonstration of some gsl reaction vectors
        
        gsl_vector *rv1 = gsl_vector_alloc (ISOTOPES);
        gsl_vector *rv2 = gsl_vector_alloc (ISOTOPES);
        gsl_vector *rv3 = gsl_vector_alloc (ISOTOPES);
        
        // Fill the vectors with values for components
        for (int i = 0; i < ISOTOPES; i++){
            gsl_vector_set (rv1, i, (i+1)*3.0);
            gsl_vector_set (rv2, i, 1.0/((i+1)*3.0));
            gsl_vector_set (rv3, i, 1.0/((i+1)*3.0));
        }
        
        // Print values of vector components
        printf("\nExample: components of GSL vector rv1:");
        for (int i = 0; i < ISOTOPES; i++){
            printf ("\nrv1_%d = %8.5f", i, gsl_vector_get (rv1, i));
        }
        printf("\n");
        
        // As test, set field values for the Species object isotope[0] and read them back
        // using the array of Species objects isotope[i].
        
        printf("\nTest of set and get functions in class Species using arrays:");
        char stg1[5] = {'1','2','C'};
        isotope[0].setisoLabel(stg1);
        isotope[0].setZ(6);
        isotope[0].setN(6);
        isotope[0].setA(12);
        isotope[0].setY(0.12);
        isotope[0].setM(0.21);
        
        printf("\n%s Z=%d N=%d A=%d Y=%g massExcess=%g\n", 
            isotope[0].getLabel(), isotope[0].getZ(), 
            isotope[0].getN(), isotope[0].getA(), 
            isotope[0].getY(), isotope[0].getM()
        );
        
        // Now set field values for the Species object isotope[1] and read them back
        // using pointers instead.  Set Species pointer to address of Species object 
        // isotope[1]
        
        SpeciesPtr = &isotope[1];
        
        printf("\nTest of set and get functions in class Species using pointers:");
        char stg2[5] = {'1','6','O'};
        
        // Use pointers to execute set() functions of Species object isotope[1]
        SpeciesPtr->setisoLabel(stg2);
        SpeciesPtr->setZ(8);
        SpeciesPtr->setN(8);
        SpeciesPtr->setA(16);
        SpeciesPtr->setY(0.13);
        SpeciesPtr->setM(0.31);
        
        // Use pointers to execute the getX() functions of Species object isotope[1]
        char tlabel[5];
        for(int k=0; k<5; k++){
            tlabel[k] = SpeciesPtr->getLabel(k);  // 2nd form of overloaded getLabel function
        }
        int Z2 = SpeciesPtr->getZ();
        int N2 = SpeciesPtr->getN();
        int A2 = SpeciesPtr->getA();
        double Y2 =  SpeciesPtr->getY();
        double M2 = SpeciesPtr->getM();
        printf("\n%s Z=%d N=%d A=%d Y=%g massExcess=%g\n", tlabel, Z2, N2, A2, Y2, M2);
        
        // As a test, set and get some fields of the Reaction object reaction [] using both
        // array notation and pointer notation.
        
        string str1 = "Now is";
        int rind = 3;
        reaction[0].setreacIndex(rind);
        int test1 = reaction[0].getreacIndex();
        cout << "\nTest of setting and getting fields of Reaction object reaction[] using arrays\n";
        cout << "reaction[0].reacIndex=" << test1 << endl;
        cout << "\nTest of setting and setting getting fields of Reaction object reaction[] with pointers\n";
        ReactionPtr = &reaction[0];
        ReactionPtr->setreacIndex(4);
        (ReactionPtr+1)->setreacIndex(5);
        cout << "Reaction index for reaction[0] is " << ReactionPtr->getreacIndex() <<
        "; reaction index for reaction[1] is " << (ReactionPtr+1)->getreacIndex() << endl;
        
        // Test of utility Utilities::returnNetIndexZN(Z,N)
        int testZ = 6;
        int testN = 6;
        int netindy = Utilities::returnNetIndexZN(testZ, testN);
        printf("\nTEST utility returnNetIndexZN(%d,%d):\n", testZ, testN);
        if(netindy < 0){
            printf("\nERROR:Species Z=%d N=%d not in network\n",testZ, testN);
        } else {
        printf("index=%d %s Z=%d N=%d\n",
            netindy,isoLabel[netindy], Z[netindy], N[netindy]);
        }
        
        // Test of Utilities::returnNetIndexSymbol(symbol)
        char testLabel[5] = "16O";
        int tester = Utilities::returnNetIndexSymbol(testLabel);
        printf("\nTEST utility returnNetIndexSymbol(%s):\n",testLabel);
        if(tester < 0){
            printf("Error: Species %s not in network\n",testLabel);
        } else {
        printf("index=%d %s Z=%d N=%d\n",
            tester,testLabel,Z[tester],N[tester]);
        }
        
        // Test of Utilities::isInNet(Z,N)
        bool isIt;
        int testz=6;
        int testn=6;
        isIt = Utilities::isInNet(testz,testn);
        printf("\nTEST Utilities::isInNet(Z,N):\nZ=%d N=%d %s",
            testz,testn,isIt ? "true" : "false");
        printf("\n");
        
        // Test of Utilities::minimumOf(i, j) for two integers
        int ti=8;
        int tj=4;
        printf("\nTEST Utilities::minimumOf(int, int):");
        printf("\n%d is the minimum of %d and %d\n", Utilities::minimumOf(ti, tj), ti, tj);
        
        // Test of utility Utilities::maximumOf (x, y) for two doubles
        double tii = -2.1;
        double tjj = -3.4;
        printf("\nTEST Utilities::maximumOf(double, double):");
        printf("\n%g is the maximum of %g and %g\n\n", 
            Utilities::maximumOf(tii, tjj), tii, tjj);
        
        // Test of Utilities::stringToChar(string)
        printf("TEST Utilities::stringToChar(string) to convert string to Char array:\n");
        string ss = "Now is the time";
        printf("Char array = %s\n\n", Utilities::stringToChar(ss));
        
        // Test of static method ReactionVector::compareGSLvectors to 
        // compare two GSL vectors. Returns 0 if not equal, 1 if equal, and
        // 2 if one vector is the negative of the other. Arguments
        // of compareGSLvectors are pointers to the two GSL vectors.
        
        int indy1 = 4;
        int indy2 = 6;
        printf("TEST compareGSLvectors() in class ReactionVector\n");
        int check = ReactionVector::compareGSLvectors(rvPt+indy1, rvPt+indy2);
        if(check==0){
            printf("Reaction vectors rv[%d] and rv[%d] are not equal: check=%d\n", 
                indy1, indy2, check);
        } else if (check==1){
            printf("Reaction vectors rv[%d] and rv[%d] are equal: check=%d\n", 
                indy1, indy2, check); 
        } else if (check==2){
            printf("Reaction vectors rv[%d] and -rv[%d] are equal: check=%d\n", 
                indy1, indy2, check);
        } 
        
        // Test of eulerUpdate(double FplusSum, double FminusSum, double Y, double dt)
        // and asymptoticUpdate(double Fplus, double Fminus, double Y, double dt)
        
        // Forward Euler
        printf("\nTEST of forward Euler updater\n");
        double fplussum = 801.3;
        double fminussum = 800.0;
        double yy = 0.22;
        double dtt = 0.0001;
        double Yupdate = Integrate::eulerUpdate(fplussum, fminussum, yy, dtt);
        printf("\nYupdate = %9.5e\n", Yupdate);
        
        // Asymptotic
        printf("\nTEST of asymptotic updater\n");
        double fplus = 901.3;
        double fminus = 900.0;
        yy = 0.21;
        dtt = 0.001;
        Yupdate = Integrate:: asymptoticUpdate(fplus, fminus, yy, dtt);
        printf("\nYupdate = %9.5e\n", Yupdate);
        
        // Test of checkAsy(double Fminus, double Y, double dt)
        printf("\nTEST of check for asymptotic condition\n");
        fminus = 2.01e6;
        yy = 0.20;
        dtt = 1.0e-7;
        bool btest = Integrate::checkAsy(fminus, yy, dtt);
        if(btest){
            printf("\nIsotope is asymptotic since ck>1\n");
        } else {
            printf("\nIsotope is NOT asymptotic since ck<1\n");
        }
        
        // Test of CPU timer by executing a long, pointless loop
        Utilities::testTimerCPU();
    
    }        // End display of test functions if showFunctionTests==1

   
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



/* Function readNetwork to read the network data file line by line, with the filename as argument.
 * This file is expected to have 4 lines per isotope with the line structure
 *	 isotopeSymbol A  Z  N  Y  MassExcess
 *	 pf00 pf01 pf02 pf03 pf04 pf05 pf06 pf07
 *	 pf10 pf11 pf12 pf13 pf14 pf15 pf16 pf17
 *	 pf20 pf21 pf22 pf23 pf24 pf25 pf26 pf27
 * where isotopeSymbol is an isotope label, A=Z+N is the atomic mass number, Z is the proton number, 
 * N is the neutron number, Y is the current abundance, MassExcess is the mass
 * excess in MeV, and the pf are 24 values of the partition function for that isotope at
 * different values of the temperature that will form a table for interpolation in temperature.
 * The assumed 24 values of the temperature for the partition function table are in array Tpf:
 * { 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 
 * 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 } in units of 10^9 K.
 * All fields on a line are separated by a blank space and there is no whitespace in the isotopeSymbol.
 * The type signature of these four lines corresponding to a single isotope is
 *	string int int int double double
 *	double double double double double double double double
 *	double double double double double double double double
 *	double double double double double double double double
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

void readNetwork (char *fileName){
    char line[60];
    char isoSymbol[5];
    int z, n, a;
    double y, mass;
    double pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7;
    
    // Open a file for reading  
    fr = fopen (fileName, "r");
    
    // Exit if the file doesn't exist or can't be read
    if( fr == NULL ){
        printf ("*** File Input Error: No readable file named %s\n",fileName);
        exit(1) ;
    }
    
    // Read in the file line by line
    
    int isoIndex = -1;
    int isoSubIndex = 3;
    
    if(displayInput==1) printf("\n--- Read in network and partition function ---\n");
    
    // Read lines until NULL encountered. Lines can contain up to 60 characters
    
    while(fgets(line, 60, fr) != NULL){
        isoSubIndex ++;
        if(isoSubIndex == 4){
            isoSubIndex = 0;
            isoIndex ++;
            
            // Read 1st line
            sscanf (line, "%s %d %d %d %lf %lf", isoSymbol, &a, &z, &n, &y, &mass);
            
            if(displayInput == 1){
                printf("\n%s %d %d %d %f %f\n", isoSymbol, a, z, n, y, mass);
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
            
            // Scan and parse a partition function line. 
            sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf", &pf0, &pf1, &pf2, &pf3, 
                    &pf4, &pf5, &pf6, &pf7);
            
            if(displayInput == 1){
                printf("%f %f %f %f %f %f %f %f\n", pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7);
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
    
}    // End of function readNetwork (char *fileName)



/* Function to read rate parameter data file line by line, with filename as argument.
 * This file is expected to have one reaction per line with the line structure
 *    p0 p1 p2 p3 p4 p5 p6 reactionLabel
 * where the pn are the values of the 7 Reaclib parameters for a reaction,
 * reactionLabel is a label for the reaction that must contain no whitespace, and
 * all fields on a line are separated by a blank space.
 */

void readLibraryParams (char *fileName)
{
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
     *     double double double double double double double string
     * each separated by a space, with no whitespace in the string.
     * (See http://stackoverflow.com/questions/2854488/reading-a-string-with-spaces-with-sscanf
     * for how to read string with spaces.)
     */
    
    int n = -1;
    int subindex = -1;
    
    if(displayInput == 1) printf("\n--- READ IN REACTIONS DATA---");
    
    // Read lines until NULL encountered. Lines can contain up to 120 characters.  In the
    // data file each reaction has 8 lines of entries.  The counter n holds the reaction number
    // and the counter subindex holds the line number for the current reaction.  The 
    // switch(subindex) determines which line we are reading (0-7) for a given reaction labeled by n.
    
    while(fgets(line, 120, fr) != NULL) {
        subindex ++;
        switch(subindex){
            
            case 0:   // 1st line
                
                n++;
                sscanf (line, "%s %d %d %d %d %d %d %d %lf %lf %d", rlabel, &i0, &i1, &i2, &i3, &i4, &i5, &i6, 
                        &sf, &q, &i7);
                
                // Store in the Reaction class instance reaction[n]
                ReactionPtr = &reaction[n];
                ReactionPtr->setreacIndex(n);
                ReactionPtr->setreacString(rlabel);   // Reaction setter will also fill reacLabel in main
                ReactionPtr->setreacGroupClass(i0);   // Reaction setter will also fill RGclass[] in main
                ReactionPtr->setRGmemberIndex(i1);    // Reaction setter will also fill RGMemberIndex[] in main
                ReactionPtr->setreacClass(i2);
                ReactionPtr->setnumberReactants(i3);  // Reaction setter will also fill NumReactingSpecies[] in main
                ReactionPtr->setnumberProducts(i4);   // Reaction setter will also fill NumProducts[] in main
                ReactionPtr->setisEC(i5);
                ReactionPtr->setisReverse(i6);
                ReactionPtr->setQ(q);
                ReactionPtr->setprefac(sf);
                ReactionPtr->setispeforward(i7);
                //isPEforward[n] = i7;
                
                if(displayInput == 1){
                    printf("\n\nReaction %d: ",n);
                    printf("%s reaclib=%d RG:(%s) RGclass=%d RGmemberIndex=%d",
                        reaction[n].getreacChar(), 
                        reaction[n].getreacIndex(),
                        reaction[n].getreacGroupChar(),
                        reaction[n].getreacGroupClass(), 
                        reaction[n].getRGmemberIndex());
                    printf("\n--------------------------------------------------------------------------------");
                    printf("\n%s %d %d %d %d %d %d %d %f %f", 
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
                
                if(displayInput == 1) printf("\n%f %f %f %f %f %f %f", 
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
                ReactionPtr -> setreactantZ(tempZN);    // Reaction setter also changes reacZ[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        printf("\n  Reactant[%d]: Z=%d", mm, reacZ[n][mm]);
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
                ReactionPtr -> setreactantN(tempZN);   // Reaction setter also changes reacN[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        printf("\n  Reactant[%d]: N=%d", mm, reacZ[n][mm]);
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
                ReactionPtr -> setproductZ(tempZN);    // Reaction setter also changes prodZ[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        printf("\n  Product[%d]: Z=%d", mm, prodZ[n][mm]);
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
                ReactionPtr -> setproductN(tempZN);    // Reaction setter also changes prodN[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        printf("\n  Product[%d]: N=%d", mm, prodN[n][mm]);
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
                ReactionPtr -> setreactantIndex(tempZN);   // setter also changes ReactantIndex[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumReactingSpecies[n]; mm++) {
                        printf("\n  ReactantIndex[%d]: N=%d", mm, ReactantIndex[n][mm]);
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
                ReactionPtr -> setproductIndex(tempZN);    // setter also changes ProductIndex[][] in main
                
                if(displayInput == 1){
                    for(int mm=0; mm<NumProducts[n]; mm++) {
                        printf("\n  ProductIndex[%d]: N=%d", mm, ProductIndex[n][mm]);
                    }
                }
                
                subindex = -1;
                
                break;
                
        }

    }
    
    // Normally numberReactions=SIZE, but numberReactions counts the actual number of
    // reactions read in.
    
    numberReactions = n+1;
    
    printf("\n%d REACTIONS\n",numberReactions);
    
    fclose(fr);           // Close the file
    
}    // End of function readLibraryParams (char *fileName)



// Function to print out the network isotopes, mass excesses, and the entries in the 
// partition function table for each isotope.

void writeNetwork()
{
    printf("\n%d ISOTOPES IN NETWORK:\n\n",numberSpecies);
    printf("Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    for (int i=0; i<numberSpecies; i++){
        printf("%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n",  
               i, isoLabel[i], AA[i], Z[i], N[i], 
               Y[i], X[i], massExcess[i]);
    }
    
    // Print out partition function table from isotope[]
    
    printf("\n\nPARTITION FUNCTION TABLE from Species object isotope[]:\n");
    printf("\n T9 = ");
    for(int k=0; k<24; k++){
        printf("%4.2f ", isotope[0].getTpf(k));
    }
    for(int i=0; i<ISOTOPES; i++){
        printf("\n");
        printf("%-5s ",isotope[i].getLabel());
        for(int j=0; j<24; j++){
            printf("%4.2f ", isotope[i].getpf(j)); 
        }
    }
    printf("\n\n");
    
}   // End of function writeNetwork()


// Function to print out all the rates. The label can be used to distinguish cases
// if called more than once.

void  writeRates(char *label)
{
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

void setRG(int index, int RGclass, int RGindex){
    reaction[index].setreacGroupClass(RGclass);
    reaction[index].setrgindex(RGindex);
}


// Helper function that can be called from another class to set some fields
// in the ReactionGroup objects RG[].

void assignRG(){
    
    printf("\n\n\nPOPULATING RG[] OBJECT FIELDS\n");
    
    for(int i=0; i<numberRG; i++){   // Loop over RGs
        RG[i] = ReactionGroup(i);
        printf("\nRG = %d", RG[i].getRGn());
        int rgindex = -1;
        for(int j=0; j<SIZE; j++){   // Loop over members of each RG
            if(RGindex[j] == i){
                rgindex ++;          // Index for member reactions in RG
                RG[i].setmemberReactions(rgindex, j);
                RG[i].setrgclass(RGclass[j]);
                RG[i].setisForward(rgindex, isPEforward[j] );
                RG[i].setRGarrayIndex(i);
                RG[i].setniso(RGclass[j]);
                
                int ck1 = RG[i].getmemberReactions(rgindex);  //reacIndex of member reaction in RG[]
                RG[i].setnumberReactants(rgindex, reaction[ck1].getnumberReactants());
                RG[i].setnumberProducts(rgindex, reaction[ck1].getnumberProducts());
                
//                 printf("\n\n +++rgindex=%d numberReactants=%d numberProducts=%d\n", 
//                     rgindex, 
//                     RG[i].getnumberReactants(rgindex),
//                     RG[i].getnumberProducts(rgindex)
//                 );
                
                int upper1 = reaction[ck1].getnumberReactants();
                int indy;
                
                // Loop over reactant isotopes
                
                for(int k=0; k<upper1; k++){
                    indy = reaction[ck1].getreactantIndex(k);
                    RG[i].setreactantIsoIndex(k, indy);
                    RG[i].setisoindex(k, indy);
                    RG[i].setisoZ(k, Z[indy]);
                    RG[i].setisoN(k, N[indy]);
                    RG[i].setisoA(k, AA[indy]);
                    RG[i].setisolabel(k, isoLabel[RG[i].getisoindex(k)]);
//                     printf("\n@@@ Reactants: k=%d isoindex=%d %s",
//                            k, RG[i].getisoindex(k), RG[i].getisolabel(k)
//                     );
                    
                    //printf("\n\n@@@ k=%d %s", k, RG[i].getisolabel(k));
                    
//                    printf("\n@@@ Reactant k=%d indy=%d Z=%d N=%d A=%d isoindex=%d\n", 
//                           k, indy, 
//                           RG[i].getisoZ(k), 
//                           RG[i].getisoN(k),
//                           RG[i].getisoA(k),
//                           RG[i].getisoindex(k)
//                    );
                }
                
                // Loop over product isotopes

                int upper2 = reaction[ck1].getnumberProducts();
                for(int k=0; k<upper2; k++){
                    indy = reaction[ck1].getproductIndex(k);
                    RG[i].setproductIsoIndex(k, indy);
                    RG[i].setisoindex(k+upper1, indy);
                    
                    RG[i].setisoZ(k+upper1, Z[indy]);
                    RG[i].setisoN(k+upper1, N[indy]);
                    RG[i].setisoA(k+upper1, AA[indy]);
                    
                    RG[i].setisolabel(k+upper1, isoLabel[RG[i].getisoindex(k+upper1)]);
//                     printf("\n@@@ Products: k=%d isoindex=%d %s",
//                            k, RG[i].getisoindex(k+upper1), RG[i].getisolabel(k+upper1)
//                     );
                    
//                     printf("\n@@@ Product k=%d indy=%d Z=%d N=%d A=%d isoindex=%d\n", 
//                            k+RG[i].getnumberReactants(rgindex), indy, 
//                            RG[i].getisoZ(k+upper1), 
//                            RG[i].getisoN(k+upper1),
//                            RG[i].getisoA(k+upper1),
//                            RG[i].getisoindex(k+upper1)
//                     );
                }
                
                RG[i].setreacString(rgindex, reaction[ck1].getreacString());
                
                printf("\nreacIndex=%d niso=%d memberIndex=%d %s RGclass=%d isForward=%d", 
                    RG[i].getmemberReactions(rgindex),
                    RG[i].getniso(),
                    rgindex,
                    RG[i].getreacString(rgindex),  
                    RG[i].getrgclass(),
                    RG[i].getisForward(rgindex)
                );
            }
            
        }
        
        RG[i].setnumberMemberReactions(rgindex+1);
        printf("\nMember reactions = %d\n", RG[i].getnumberMemberReactions());
    }
}       // End function assignRG()


// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesfplus(int index, double fp){
    isotope[index].setfplus(fp);
}

// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesfminus(int index, double fp){
    isotope[index].setfminus(fp);
}

// Helper function that can be called from another class to set field
// in the Species objects isotope[].

void setSpeciesdYdt(int index, double dydt){
    isotope[index].setdYdt(dydt);
}
