/*
 * Code to implement explicit algebraic integration of astrophysical thermonuclear networks.
 * Execution assuming use of Fedora Linux: Compile with
 *     gcc explicitMatrix.c -o explicitMatrix -lgsl -lgslcblas -lm -lgslcblas -lm -lstdc++
 * Resulting compiled code can be executed with
 *     ./explicitMatrix  | tee temp.txt
 * where | tee temp.txt is unix shell script outputting to screen and also piped to a file temp.txt. 
 * Execution for other Linux systems, or Mac or PC, will depend on the C/C++ compiler installed on 
 * your machine but should be similar.  
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
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

// NOTE:  For different networks you must change the values of SIZE, ISOTOPES, rateLibraryFile[],
//        and networkFile[] to appropriate values.

// SIZE defines the number of reactions to be calculated.  Need min of SIZE=4395 for 365-element network, 
// 1603 for 150-isotope network, 2234 for the 194-isotope network, 3176 for the 268-isotope network,
// 1566 for the nova134 network, 48 for the alpha network, 8 for 3-alpha 4He-12C-16O network, and 28 for 
// the 7-isotope pp network. ISOTOPES defines the number of isotopes in each network (for example, ISOTOPES=7
// for the 7-isotope pp-chain network.  These sizes are hardwired for now but eventually we want to read
// them in and assign them dynamically.

#define SIZE 48                        // Max number of reactions (48 for alpha network)
#define ISOTOPES 16                    // Max isotopes in network (16 for alpha network) 16

#define LABELSIZE 35                  // Max size of reaction string a+b>c in characters
#define PF 24                         // Number entries in partition function table for each isotope
#define THIRD 0.333333333333333

// Define some CPU timing utilities. Usage:
//     START_CPU;
//     ... code to be timed ...
//     STOP_CPU;
//     PRINT_CPU
// These may be used to time CPU processes.

clock_t startCPU, stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("\nTimer: %f ms used by CPU\n",1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));

// File pointer for data read-in
FILE *fr;            

// Filename for input rates library data. The file rateLibrary.data output by the Java code through 
// the stream toRateData has the expected format for this file.  Standard test cases: 
// rateLibrary_alpha.data, rateLibrary_150.data, rateLibrary_365.data, rateLibrary_nova134.data,
// rateLibrary_3alpha.data, rateLibrary_pp.data.

char rateLibraryFile[] = "rateLibrary_alpha.data";  //"rateLibrary_alpha.data";

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format for 
// this file. Standard test cases: CUDAnet_alphasolar.inp, CUDAnet_150solar.inp,
// CUDAnet_365solar.inp, CUDAnet_nova134.inp, CUDAnet_3alpha.inp, CUDAnet_pp.inp.

char networkFile[] = "CUDAnet_alpha.inp"; //"CUDAnet_alpha.inp";

// Control some of the printout (1 to include, 0 to suppress)
int displayInput = 1;
int ratesCPU = 1;

// Control diagnostic printout of details (1 to print, 0 to suppress)
int showParsing = 1;
int showFparsing = 1;
int showFluxCalc = 1;

// Function Signatures:
void devcheck(int);
void readLibraryParams(char *);
void readNetwork(char *);
void writeNetwork(void);
void testTimerCPU(void);
void computeRatesCPU(double);
void computeFluxCPU(void);
void writeRates(char *);
void writeAbundances(void);
void parseF(void);
int returnNetIndexZN(int,int);
int returnNetIndexSymbol(char *);
bool isInNet(int,int);
int minimumOf(int,int);
int maximumOf(int,int);
void makeReactionVectors(void);
int compareGSLvectors(gsl_vector*, gsl_vector*);
void sortReactionGroups(void);

// Variables to hold data that will be read in

double P0[SIZE];       // Array holding library rate parameters p0
double P1[SIZE];       // Array holding library rate parameters p1
double P2[SIZE];       // Array holding library rate parameters p2
double P3[SIZE];       // Array holding library rate parameters p3
double P4[SIZE];       // Array holding library rate parameters p4
double P5[SIZE];       // Array holding library rate parameters p5
double P6[SIZE];       // Array holding library rate parameters p6
double Q[SIZE];        // Array holding library entry for reaction Q-value
double Rate[SIZE];     // Array holding the computed rate for each reaction
double Flux[SIZE];     // Array holding the computed fluxes for each reaction

char reacLabel[SIZE][LABELSIZE]; // Array holding reaction labels (e.g. he4+c12-->o16)
int reaclibClass[SIZE];          // Reaclib class index for reaction (1-8) 
int RGclass[SIZE];               // Reaction Group class (PE) for reaction (1-5)
int RGmemberIndex[SIZE];         // Member index within reaction group
int NumReactingSpecies[SIZE];    // Number of reactant isotopes for the reaction
int NumProducts[SIZE];           // Number of product isotopes for the reaction
int isEC[SIZE];                  // Whether electron capture (0=false; 1=true)
int isReverseR[SIZE];            // Whether inverse reaction (matters for partition func)
double Prefac[SIZE];             // The statistical prefactor for each reaction

int reactantZ[SIZE][4];          // Holds Z for each reactant isotope
int reactantN[SIZE][4];          // Holds N for each reactant isotope
int productZ[SIZE][4];           // Holds Z for each product isotope
int productN[SIZE][4];           // Holds N for each product isotope
int ReactantIndex[SIZE][4];      // Index of isotope vector for each reactant isotope
int ProductIndex[SIZE][4];       // Index of isotope vector for each product isotope
int Reactant1[SIZE];             // Isotope index of 1st reactant
int Reactant2[SIZE];             // Isotope index of 2nd reactant (if 2-body or 3-body)
int Reactant3[SIZE];		     // Isotope index of 3rd reactant (if 3-body)

// ------ Species data in the following 7 arrays also contained in struct networkData

 int Z[ISOTOPES];                 // Array holding Z values for isotopes
 int N[ISOTOPES];                 // Array holding N values for isotopes
 double AA[ISOTOPES];             // Array holding A values for isotopes
 double Y[ISOTOPES];              // Array holding abundances Y for isotopes
 double X[ISOTOPES];              // Array holding mass fractions X for isotopes
 double massExcess[ISOTOPES];     // Array holding mass excesses for isotopes
 char isoLabel[ISOTOPES][5];      // Isotope labels (max length 5 characters; e.g. 238pu)
// -----------

int numberSpecies;                // # of species in network (generally = ISOTOPES)
int numberReactions;              // # of reactions in the network (generally = SIZE)

// Temperatures in units of 10^9 K for partition function table (see pf[][]). These are
// copied from the corresponding variable Tpf[] in the Java code.

const double Tpf[] = { 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f,
    1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f };
    
// Array holding partition function values for the 24 temperatures given in Tpf[]
// for each isotope

double pf[ISOTOPES][PF]; 

// Dens[] will hold the density factor for fluxes. Dens[0] = 1.0 (1-body reactions), 
// Dens[1]=rho (2-body reactions), Dens[2]=rho*rho (3-body reactions)

double Dens[3];  

// Array with entries +1 if a reaction increases the population of the isotope (contributes to 
// F+), -1 if it decreases it (contributes to F-) and 0 if the reaction does not change the population
// of the isotope. This array is populated in the function parseF().  It is characteristic of the
// structure of the network and thus has to be calculated only once for a given network.

int reacMask[ISOTOPES][SIZE];

// Define an integer reaction vector array RV[reactionNumber][species number] where reactNumber
// labels the reaction and species number labels the components of the vector

int RV[SIZE][ISOTOPES]; 

// Define an array rv[] and corresponding pointers that will hold GSL vectors corresponding
// to the reaction vectors for the system.  This will be implemented in the function
// makeReactionVectors().

gsl_vector rv[SIZE];   // Array of type gsl_vector to hold GSL vectors
gsl_vector *rvPt;      // Pointer to rv[] array

int numberRG;          // Number of partial equilibrium reaction groups

// Define array to hold the reaction group index for each reaction.  This array is
// populated by the function sortReactionGroups()

int RGindex[SIZE];

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
// Fminus arrays. Corresponding memory will be allocated dynamically below malloc

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
 * These will be allocated dynamically below.
 */

int* MapFplus;    // Index mapper for Fplus (Dim totalFplus)
int* MapFminus;   // Index mapper for Fminus (Dim totalFminus)

// Arrays for temporary storage. Will be allocated dynamically below
int* tempInt1;
int* tempInt2;


/* Create a struct networkData to hold the info on species of the network, with the 
 number of entries stored in the struct given by the integer ISOTOPES.  Access fields 
 for each isotope with the pointer networkDataPtr like so:
 
        char labelfield[5];         // isotope label
        for(int k=0; k<5; k++){
            labelfield[k] = (networkDataPtr+i)->isoLabel[k];  
        }
        int Afield = (networkDataPtr+i)->A;              // mass number
        int Zfield = (networkDataPtr+i)->Z;              // proton number
        int Nfield = (networkDataPtr+i)->N;              // neutron number
        double Yfield = (networkDataPtr+i)->Y;           // abundance
        double Mfield = (networkDataPtr+i)->massExcess;  // mass excess
        
 where i numbers the species in the network (with i=0 for the first entry and 
 i=ISOTOPES-1 for the last entry). For example, 
 
        int Zfield = (networkDataPtr+3)->Z;
        
 would access the proton number Z of the 4th species in the network (numbered from zero) and
 assign it to the integer variable Zfield.  The order of species stored in struct is 
 determined by the order of data read in using the function readNetwork(char *file), where 
 the struct is populated as data are read in.
 */

struct networkData{
    char isoLabel[5];   // symbol for isotope
    int Z;              // proton number
    int N;              // neutron number
    int A;              // atomic mass number
    double Y;           // abundance
    double massExcess;  // mass excess
};

// Pointer used to access the struct networkData
struct networkData *networkDataPtr;


// Add a C++ class to store species data and functions.  Need to save with .cpp 
// extension and compile with gcc flag -lstdc++ as follows:
//    gcc explicitMatrix.cpp -o explicitMatrix -lgsl -lgslcblas -lm -lstdc++
// (assuming Fedora Linux)

class Species{
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
    private:
        
        char isoLabel[5];   // symbol for isotope
        int Z;              // proton number
        int N;              // neutron number
        int A;              // atomic mass number
        double Y;           // abundance
        double massExcess;  // mass excess
    
    public:
        
        // Public setter function to set values of private class fields
        
        void set (char *lpt, int z, int n, int a, double y, double m){
            
            // Fill character array.  *lpt points to initial character.
            // All characters of the character array are accessed by
            // incrementing *lpt.
            
            for(int j=0; j<5; j++){
                isoLabel[j] = *(lpt+j);
            }
            Z = z;
            N = n;
            A = a;
            Y = y;
            massExcess = m;
        };
        
        // Getter functions to return values of class fields
        
        int getZ() {return Z; };
        int getN() {return N; };
        int getA() {return A; };
        double getY() {return Y; };
        double getM(){return massExcess; };
        
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
        
        char *getLabel() {return isoLabel; };         // return pointer to label array
        char getLabel(int k) {return isoLabel[k]; };  // return kth character of array
};


// Create an array of Species objects isotope[] to hold information and functions
// for the isotopic species in the network. Each element of the array will be
// a Species object corresponding to a different isotope of the network.

Species isotope[ISOTOPES];



// Main CPU routine

int main()
{  
    // Allocate memory for struct networkData with pointer networkDataPtr
    networkDataPtr = (struct networkData*) malloc(ISOTOPES * sizeof(struct networkData));
    
    // As test, set field values for the Species object isotope[0] and read them back
    // using the array isotope[i]
    
    printf("\nTest of set and get functions in class Species using arrays:");
    char stg1[5] = {'1','2','C'};
    isotope[0].set (stg1, 6, 6, 12, 0.12, 0.21);
    printf("\n%s Z=%d N=%d A=%d Y=%g massExcess=%g\n", isotope[0].getLabel(), isotope[0].getZ(), 
        isotope[0].getN(), isotope[0].getA(), isotope[0].getY(), isotope[0].getM()
    );
    
    // Now set field values for the Species object isotope[1] and read them back
    // using pointers
    
    // Declare pointer used to access the fields and functions of class Species
    Species *SpeciesPtr;
    
    // Set Species pointer to address of Species object isotope[1]
    SpeciesPtr = &isotope[1];
    
    printf("\nTest of set and get functions in class Species using pointers:");
    char stg2[5] = {'1','6','O'};
    
    // Use pointers to execute the set() function of Species object isotope[1]
    SpeciesPtr->set(stg2, 8, 8, 16, 0.13, 0.31);
    
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
    
    
    // Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
    // realistic calculation the temperature and density will be passed from the hydro 
    // code in an operator-split coupling of this network to hydro. Here we hardwire
    // them for testing purposes.  These will be used to calculate the reaction
    // rates in the network. Since we are assuming operator splitting, the temperature 
    // and density are assumed constant for each network integration.
    
    double T9 = 6.0f;
    double rho = 1.0e8;
    
    // Set the range of time integration and the initial timestep.  In an operator-split
    // coupling tmax will come from the hydro and dt_init will likely be the last timestep
    // of the previous network integration (for the preceding hydro timestep). Here we
    // hardwire them for testing purposes.
    
    double tmax = 1e-11;
    double dt_init = 1e-17; 
    
    // Demonstration of some gsl reaction vectors
    
    gsl_vector *rv1 = gsl_vector_alloc (ISOTOPES);
    gsl_vector *rv2 = gsl_vector_alloc (ISOTOPES);
    gsl_vector *rv3 = gsl_vector_alloc (ISOTOPES);
    
    // Fill the vectors with values for components
    for (int i = 0; i < ISOTOPES; i++)
    {
        gsl_vector_set (rv1, i, (i+1)*3.0);
        gsl_vector_set (rv2, i, 1.0/((i+1)*3.0));
        gsl_vector_set (rv3, i, 1.0/((i+1)*3.0));
    }
    
    // Print values of vector components
    printf("\nExample: components of GSL vector rv1:");
    for (int i = 0; i < ISOTOPES; i++)
    {
        printf ("\nrv1_%d = %8.5f", i, gsl_vector_get (rv1, i));
    }
    printf("\n");
    
    // Read in network file and associated partition functions.  This is required only
    // once at the beginning of the entire calculation.  
    
    char *networkFilePtr = networkFile;
    readNetwork(networkFilePtr);
    writeNetwork();

    // Write out fields of struct networkData
    printf("\n\nFIELDS OF STRUCT networkData{isoLabel, Z, N, A, Y, massExcess}\n");
    printf("that stores info on each species in the network:\n");
    for(int i=0; i<ISOTOPES; i++){
        
        // The struct networkData now holds the data on isotopic species read in
        // from the data file.  Use pointers to retrieve data from struct for ith species
        
        char labelfield[5];         // istope label
        for(int k=0; k<5; k++){
            labelfield[k] = (networkDataPtr+i)->isoLabel[k];  
        }
        int Zfield = (networkDataPtr+i)->Z;              // proton number
        int Nfield = (networkDataPtr+i)->N;              // neutron number
        int Afield = (networkDataPtr+i)->A;              // mass number
        double Yfield = (networkDataPtr+i)->Y;           // abundance
        double Mfield = (networkDataPtr+i)->massExcess;  // mass excess
        
        printf("\ni=%d  isoLabel=%s  Z=%d  N=%d  A=%d  Y=%6.4f  MassExcess=%5.2lf", 
            i,
            labelfield,
            Zfield,
            Nfield,
            Afield,
            Yfield,
            Mfield 
        );
    }
    printf("\n\nAccess field of species i (counting from 0) with pointer: (networkDataPtr+i)->field");
    printf("\nEXAMPLE: int myZ = (networkDataPtr+0)->Z to assign proton number Z of 1st (i=0) species in network to variable myZ\n");
    printf("\n");
    
    printf("\n\nDATA FIELDS of Species object isotope[] after network data read in\n");
    for(int i=0; i<ISOTOPES; i++){
        printf("\nisotope[%2d]: %5s  Z=%2d  N=%2d  A=%3d  Y=%6.4f  massExcess=%6.4f", 
               i, isotope[i].getLabel(), isotope[i].getZ(), isotope[i].getN(), 
               isotope[i].getA(), isotope[i].getY(), isotope[i].getM());
    }
    printf("\n\n");
    
    
    // Read in rate library data from a file. This is required only once, at the
    // beginning of the entire calculation.
    
    char *rateLibraryFilePtr = rateLibraryFile;
    readLibraryParams(rateLibraryFilePtr);
    
    // Multiply the prefactor by the appropriate density factors (1 for 1-body,
    // rho for 2-body, and rho^2 for 3-body. This is required at the beginning of
    // each network integration of the hydro timestep, since the density will generally
    // change over a hydro timestep in each zone.
    
    Dens[0] = 1.0f;
    Dens[1] = rho;
    Dens[2] = rho*rho;
    
    for(int i=0; i<SIZE; i++)
    {
        Prefac[i] *= (Dens[NumReactingSpecies[i]-1]);
    }
    
    // Compute the temperature-dependent ReacLib rates 
    if(ratesCPU == 1){
        
//      // Test the CPU timer by executing a long, pointless loop
//      testTimerCPU();
         
        START_CPU     // Start a timer for rate calculation
        
        // Compute rates for all reactions for this temperature
        computeRatesCPU(T9);
        
        STOP_CPU;     // Stop the timer
        
        // Display the rates. (Note the two different techniques used in calling
        // writeRates here and for "GPU" in original CUDA code). Here we use a 
        // pointer; in the original CUDA kernel we passed an array.)
        
        char label[] = "on CPU";
        char *labelPtr = label;
        writeRates(labelPtr);
        
        PRINT_CPU;    // Print timing information for rate calculation
    }
    
    // Compute total flux for each reaction
    
    printf("\n--- TOTAL FLUXES (T9=%4.2f, rho=%7.3e) ---\n",T9,rho);
    
    computeFluxCPU();
    
    // Find for each isotope all reactions that change its population.  This analysis of
    // the network is required only once at the very beginning of the calculation (provided
    // that the network species and reactions remain the same for the entire calculation).
    // The work is done by the function parseF().
    
    // Number of F+ and F- components for each isotope
    numFluxPlus = (int*) malloc(sizeof(int) * numberSpecies);
    numFluxMinus = (int*) malloc(sizeof(int) * numberSpecies);
    
    // Arrays for temporary storage
    tempInt1 = (int*) malloc(sizeof(int) * numberSpecies * numberReactions/2);
    tempInt2 = (int*) malloc(sizeof(int) * numberSpecies * numberReactions/2);
    
    /*
     * Use parseF() to find the contributions to F+ and F- of each reaction for each isotope.  
     * This is executed only once at the beginning of the entire calculation to determine 
     * the structure of the network.
     */
    
    parseF();
    
    // Print out the network species vector
    printf("\n\nNETWORK SPECIES VECTOR (%d components):\n\nIndex  Species    Z     N",numberSpecies);
    for(int i=0; i<numberSpecies; i++){
        printf("\n%5d    %5s  %3d  %4d",i,isoLabel[i],Z[i], N[i]);
    }
    printf("\n");
    
    // Use the information gleaned from parseF() to define the reaction vectors
    // for the network.
    
    makeReactionVectors();
    
    // Allocate 1D arrays to hold non-zero F+ and F- for all reactions for all isotopes,
    // the arrays holding the species factors FplusFac and FminusFac, 
    // and also arrays to hold their sums for each isotope. Note that parseF() must
    // be run first because it determines totalFplus and totalFminus.
    
    Fplus = (double*) malloc(sizeof(double) * totalFplus);
    Fminus = (double*) malloc(sizeof(double) * totalFminus);
    FplusFac = (double*) malloc(sizeof(double) *totalFplus);
    FminusFac = (double*) malloc(sizeof(double) * totalFminus);
    FplusSum = (double*) malloc(sizeof(double) * numberSpecies);
    FminusSum = (double*) malloc(sizeof(double) * numberSpecies);
    
    // Arrays that hold the index of the boundary between different isotopes in the
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
    
    FplusIsotopeCut[0] = numFluxPlus[0];
    FminusIsotopeCut[0] = numFluxMinus[0];
    for(int i=1; i<numberSpecies; i++)
    {
        FplusIsotopeCut[i] = numFluxPlus[i] + FplusIsotopeCut[i-1];
        FminusIsotopeCut[i] = numFluxMinus[i] + FminusIsotopeCut[i-1];
    }
    
    int currentIso = 0;
    for(int i=0; i<totalFplus; i++)
    {
        FplusIsotopeIndex[i] = currentIso;
        if(i == (FplusIsotopeCut[currentIso]-1)) currentIso ++;
    }
    
    currentIso = 0;
    for(int i=0; i<totalFminus; i++)
    {
        FminusIsotopeIndex[i] = currentIso;
        if(i == (FminusIsotopeCut[currentIso]-1)) currentIso ++;
    }

    // Optional diagnostic output if showFparsing==1
    if(showFparsing == 1)
    {
        printf("\n\n--- MAX F+ and F- INDEX FOR EACH ISOTOPE ---\n");	
        for(int i=0; i<numberSpecies; i++)
        {
            printf("\nIsotope index = %d  %s  Max index F+ = %d  Max index F- = %d", 
                   i, isoLabel[i], FplusIsotopeCut[i]-1, FminusIsotopeCut[i]-1);
        }
    }
    
    for(int i=0; i<totalFplus; i++)
    {
        MapFplus[i] = tempInt1[i];
    }
    
    for(int i=0; i<totalFminus; i++)
    {
        MapFminus[i] = tempInt2[i];
    }
       
    // Populate the FplusMin and FplusMax arrays
    FplusMin[0] = 0;
    FplusMax[0] = numFluxPlus[0]-1;

    for(int i=1; i<numberSpecies; i++)
    {
        FplusMin[i] = FplusMax[i-1] + 1;
        FplusMax[i] = FplusMin[i] + numFluxPlus[i] -1 ;	
    }

    // Populate the FminusMin and FminusMax arrays
    FminusMin[0] = 0;
    FminusMax[0] = numFluxMinus[0]-1;
    for(int i=1; i<numberSpecies; i++)
    {
        FminusMin[i] = FminusMax[i-1] + 1;
        FminusMax[i] = FminusMin[i] + numFluxMinus[i] -1 ;	
    }
    
    // Populate the FplusFac and FminusFac arrays that hold the factors counting the
    // number of occurrences of the species in the reaction.  Note that this can only
    // be done after parseF() has been run to give reacMask[i][j].
    
    int tempCountPlus = 0;
    int tempCountMinus = 0;
    for(int i=0; i<ISOTOPES; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(reacMask[i][j] > 0)
            {
                FplusFac[tempCountPlus] = (double)reacMask[i][j];
// 				printf("\n F+  tempCountPlus=%d i=%d j=%d FplusFac=%3.1f", tempCountPlus, 
// 					   i, j, FplusFac[tempCountPlus]);
                tempCountPlus ++;
            }
            else if(reacMask[i][j] < 0)
            {
                FminusFac[tempCountMinus] = -(double) reacMask[i][j];
// 				printf("\n F-  tempCountMinus=%d i=%d j=%d FminusFac=%3.1f", tempCountMinus, 
// 					   i, j, FminusFac[tempCountMinus]);
                tempCountMinus ++;
            }	
        }
    }
    
    // Optional diagnostic output
    
    if(showFparsing == 1)
    {
        printf("\n\n\n--- %d NON-VANISHING F+ SOURCE TERMS ---\n", totalFplus);
        printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
            isoLabel[FplusIsotopeIndex[0]], FplusIsotopeIndex[0],
            numFluxPlus[FplusIsotopeIndex[0]]);
        for(int i=0; i<totalFplus; i++)
        {
            printf("\n   Isotope index = %d F+ index = %d Reac index = %d  %s", 
                   FplusIsotopeIndex[i], i, MapFplus[i], reacLabel[MapFplus[i]]); 
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
               numFluxMinus[FminusIsotopeIndex[0]]
        );
        for(int i=0; i<totalFminus; i++)
        {
            printf("\n   Isotope index = %d F- index = %d Reac index=%d  %s", 
                   FminusIsotopeIndex[i], i, MapFminus[i], reacLabel[MapFminus[i]]);
            if(i == (FminusIsotopeCut[FminusIsotopeIndex[i]] - 1) && i != totalFminus-1 )
            {
                printf("\n");
                printf("\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):", 
                       isoLabel[FminusIsotopeIndex[i+1]], FminusIsotopeIndex[i+1],
                       numFluxMinus[FminusIsotopeIndex[i+1]]
                );
            }
        }
        printf("\n\n");
    }
    
    // Populate the F+ and F- arrays from the master Flux array
    
    printf("\n--- Initial values of F+ ---\n\n");
    for(int i=0; i<totalFplus; i++)
    {
        int indy = MapFplus[i];
        Fplus[i] = FplusFac[i]*Flux[indy];
        printf("i=%d FplusFac=%3.1f Flux=%7.3e Fplus=%7.3e\n", i,FplusFac[i],Flux[indy],Fplus[i]);
    }
    
    printf("\n\n--- Initial values of F- ---\n\n");
    
    for(int i=0; i<totalFminus; i++)
    {
        int indy = MapFminus[i];
        Fminus[i] = FminusFac[i]*Flux[indy];
        printf("i=%d FminusFac=%3.1f Flux=%7.3e Fminus=%7.3e\n", i,FminusFac[i],Flux[indy],Fminus[i]);
    }
    
//     // Print out the network species vector
//     printf("\nNetwork species vector (%d components):\n\nIndex Species  Z   N",numberSpecies);
//     for(int i=0; i<numberSpecies; i++){
//         printf("\n%d     %s      %d  %d",i,isoLabel[i],Z[i], N[i]);
//     }
//     printf("\n");
    
    // Test of utility returnNetIndexZN(Z,N)
    int testZ = 8;
    int testN = 8;
    int netindy = returnNetIndexZN(testZ, testN);
    printf("\n\nTEST utility returnNetIndexZN(%d,%d):\n", testZ, testN);
    if(netindy < 0){
        printf("Error:Species Z=%d N=%d not in network\n",testZ, testN);
    } else {
    printf("index=%d %s Z=%d N=%d\n",
           netindy,isoLabel[netindy], Z[netindy], N[netindy]);
    }
    
    // Test of utility returnNetIndexSymbol(symbol)
    char testLabel[5] = "12C";
    int tester = returnNetIndexSymbol(testLabel);
    printf("\nTEST utility returnNetIndexSymbol(%s):\n",testLabel);
    if(tester < 0){
        printf("Error: Species %s not in network\n",testLabel);
    } else {
    printf("index=%d %s Z=%d N=%d\n",
           tester,testLabel,Z[tester],N[tester]);
    }
    
    // Test of utility isInNet(Z,N)
    bool isIt;
    int testz=6;
    int testn=6;
    isIt=isInNet(testz,testn);
    printf("\nTEST utility isInNet(Z,N):\nZ=%d N=%d %s\n",
           testz,testn,isIt ? "true" : "false");
    printf("\n");
    
    // Test of utility compareGSLvectors to compare two
    // GSL vectors. Returns 0 if not equal, 1 if equal, and
    // 2 if one vector is the negative of the other. Arguments
    // of compareGSLvectors are pointers to the two GSL vectors.
    
    int indy1 = 4;
    int indy2 = 6;
    
    printf("TEST utility compareGSLvectors()\n");
    int check = compareGSLvectors(rvPt+indy1, rvPt+indy2);
    if(check==0){
        printf("Reaction vectors rv[%d] and rv[%d] are not equal: check=%d\n", indy1, indy2, check);
    } else if (check==1){
        printf("Reaction vectors rv[%d] and rv[%d] are equal: check=%d\n", indy1, indy2, check); 
    } else if (check==2){
        printf("Reaction vectors rv[%d] and -rv[%d] are equal: check=%d\n", indy1, indy2, check);
    } 
    printf("\n");
    
    // Test of sortReactionGroups() to sort reactions into reaction groups by
    // comparing reaction vectors
    
    sortReactionGroups();
    printf("\n");
    
    
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
    free(networkDataPtr);  // Is this right way to free struct memory?
    
}  // End of main routine




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
    char reaction[LABELSIZE];
    double p0, p1, p2, p3, p4, p5, p6, q, sf;
    int i0, i1, i2, i3, i4, i5, i6;
    int ii[6];
    
    // Open a file for reading  
    fr = fopen (fileName, "r");
    
    // Exit if the file doesn't exist or can't be read
    if( fr == NULL )
    {
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
    
    // Read in lines until NULL encountered. Lines can contain up to 120 characters
    
    while(fgets(line, 120, fr) != NULL)
    {
        subindex ++;
        switch(subindex){
            
            case 0:
                n++;
                sscanf (line, "%s %d %d %d %d %d %d %d %lf %lf", reaction, &i0, &i1, &i2, &i3, &i4, &i5, &i6, 
                        &sf, &q);
                for(int j=0; j<LABELSIZE; j++){
                    reacLabel[n][j] = reaction[j];
                }
                
                RGclass[n] = i0;
                RGmemberIndex[n] = i1;
                reaclibClass[n] = i2;
                NumReactingSpecies[n] = i3;
                NumProducts[n] = i4;
                isEC[n] = i5;
                isReverseR[n] = i6;
                Prefac[n] = sf;
                Q[n] = q;
                
                if(displayInput == 1) printf("\n\nReaction Index = %d",n);
                if(displayInput == 1) printf("\nisReverseR = %d reaclibIndex = %d",isReverseR[n],reaclibClass[n]);
                if(displayInput == 1) printf("\n%s %d %d %d %d %d %d %d %f %f", 
                    reacLabel[n], 
                    RGclass[n],
                    RGmemberIndex[n],
                    reaclibClass[n],
                    NumReactingSpecies[n],
                    NumProducts[n],
                    isEC[n],
                    isReverseR[n],
                    Prefac[n],
                    Q[n]
                );
                
                break;
                
            case 1:
                sscanf (line, "%lf %lf %lf %lf %lf %lf %lf", &p0, &p1, &p2, &p3, &p4, &p5, &p6);
                P0[n] = p0;
                P1[n] = p1;
                P2[n] = p2;
                P3[n] = p3;
                P4[n] = p4;
                P5[n] = p5;
                P6[n] = p6;
                
                if(displayInput == 1) printf("\n%f %f %f %f %f %f %f", 
                    P0[n], 
                    P1[n],
                    P2[n],
                    P3[n],
                    P4[n],
                    P5[n],
                    P6[n]
                );
                
                break;
                
            case 2:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumReactingSpecies[n]; mm++)
                {
                    reactantZ[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  Reactant[%d]: Z=%d", mm, reactantZ[n][mm]);
                }
                
                break;
                
            case 3:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumReactingSpecies[n]; mm++)
                {
                    reactantN[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  Reactant[%d]: N=%d", mm, reactantN[n][mm]);
                }
                
                break;
                
            case 4:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumProducts[n]; mm++)
                {
                    productZ[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  Product[%d]: Z=%d", mm, productZ[n][mm]);
                }
                
                break;
                
            case 5:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumProducts[n]; mm++)
                {
                    productN[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  Product[%d]: N=%d", mm, productN[n][mm]);
                }
                
                break;
                
            case 6:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumReactingSpecies[n]; mm++)
                {
                    ReactantIndex[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  ReactantIndex[%d]: N=%d", mm, ReactantIndex[n][mm]);
                }
                
                break;
                
            case 7:
                sscanf (line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
                for(int mm=0; mm<NumProducts[n]; mm++)
                {
                    ProductIndex[n][mm] = ii[mm];
                    if(displayInput == 1) printf("\n  ProductIndex[%d]: N=%d", mm, ProductIndex[n][mm]);
                }

                subindex = -1;
                
                break;
                
        }
        
    }
    numberReactions = n+1;
    printf("\n\nnumberReactions=%d",numberReactions);
    
    for(int i=0; i<numberReactions; i++)
    {
        Reactant1[i] = ReactantIndex[i][0];
        Reactant2[i] = ReactantIndex[i][1];
        Reactant3[i] = ReactantIndex[i][2];
    }
    
    fclose(fr);           // Close the file
    
}    // End of function readLibraryParams (char *fileName)


/* Function to read the network data file line by line, with the filename as argument.
 * This file is expected to have 4 lines per isotope with the line structure
 *	 isotopeSymbol A  Z  N  Y  MassExcess
 *	 pf00 pf01 pf02 pf03 pf04 pf05 pf06 pf07
 *	 pf10 pf11 pf12 pf13 pf14 pf15 pf16 pf17
 *	 pf20 pf21 pf22 pf23 pf24 pf25 pf26 pf27
 * where isotopeSymbol is an isotope label, A=Z+N is the atomic mass number, Z is the proton number, 
 * N is the neutron number, Y is the current abundance, MassExcess is the mass
 * excess in MeV, and the pf are 24 values of the partition function for that isotope at
 * different values of the temperature that will form a table for interpolation in temperature.
 * The assumed 24 values of the temperature for the partition function table are
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

void readNetwork (char *fileName)
{
    char line[60];
    char isoSymbol[5];
    int z, n, a;
    double y, mass;
    double pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7;
    
    // Open a file for reading  
    fr = fopen (fileName, "r");
    
    // Exit if the file doesn't exist or can't be read
    if( fr == NULL )
    {
        printf ("*** File Input Error: No readable file named %s\n",fileName);
        exit(1) ;
    }
    
    // Read in the file line by line
    
    int isoIndex = -1;
    int isoSubIndex = 3;
    
    if(displayInput==1) printf("\n--- Read in network and partition function ---\n");
    
    // Read in lines until NULL encountered. Lines can contain up to 60 characters
    
    while(fgets(line, 60, fr) != NULL)
    {
        isoSubIndex ++;
        if(isoSubIndex == 4){
            isoSubIndex = 0;
            isoIndex ++;
            // Scan and parse a title line using pointers to write into fields of the struct
            // networkData
            sscanf (line, "%s %d %d %d %lf %lf", 
                (networkDataPtr+isoIndex)->isoLabel, 
                &(networkDataPtr+isoIndex)->A, 
                &(networkDataPtr+isoIndex)->Z, 
                &(networkDataPtr+isoIndex)->N, 
                &(networkDataPtr+isoIndex)->Y, 
                &(networkDataPtr+isoIndex)->massExcess);

                // printf("\n\n*** struct networkData: isoSymbol=%s A=%d Z=%d N=%d Y=%6.4lf MassExcess=%lf\n", 
                //        (networkDataPtr+isoIndex)->isoLabel, 
                //        (networkDataPtr+isoIndex)->A,
                //        (networkDataPtr+isoIndex)->Z,
                //        (networkDataPtr+isoIndex)->N,
                //        (networkDataPtr+isoIndex)->Y,
                //        (networkDataPtr+isoIndex)->massExcess
                // );
            
            // Old version for scanning directly into arrays instead of struct networkData
            //sscanf (line, "%s %d %d %d %lf %lf", isoSymbol, &a, &z, &n, &y, &mass);

            // Store variables in arrays by copying from struct networkData
            for(int k=0; k<5; k++){
                isoSymbol[k] = (networkDataPtr+isoIndex)->isoLabel[k];  
            }
            a = (int)(networkDataPtr+isoIndex)->A;
            z = (int)(networkDataPtr+isoIndex)->Z;
            n = (int)(networkDataPtr+isoIndex)->N;
            y = (double)(networkDataPtr+isoIndex)->Y;
            mass = (double)(networkDataPtr+isoIndex)->massExcess;
            if(displayInput == 1)
            {
                printf("\n%s %d %d %d %lf %lf\n", isoSymbol, a, z, n, y, mass);
            }
            Z[isoIndex] = z;
            N[isoIndex] = n;
            AA[isoIndex] = (double)a;
            Y[isoIndex] = y;
            X[isoIndex] = AA[isoIndex]*Y[isoIndex];
            massExcess[isoIndex] = mass;
            for(int j=0; j<5; j++){
                isoLabel[isoIndex][j] = (networkDataPtr+isoIndex)->isoLabel[j];
            }
            
            // Also copy data to the Species object isotope[]
            // The first argument is a pointer to the beginning of the label field
            // in the struct networkData
            
            isotope[isoIndex].set((networkDataPtr+isoIndex)->isoLabel,z,n,a,y,mass);
            //isotope[isoIndex].set(z,n,a,y,mass);
            //isotope[isoIndex].setLabel((networkDataPtr+isoIndex)->isoLabel);
            
        } else {
            
            
            
            // Scan and parse a partition function line. 
            sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf", &pf0, &pf1, &pf2, &pf3, 
                &pf4, &pf5, &pf6, &pf7);
            if(displayInput == 1)
            {
                printf("%f %f %f %f %f %f %f %f\n", pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7);
            }
            // Store the partition function table values
            int tin = isoSubIndex-1;
            pf[isoIndex][8*(tin)] = pf0;
            pf[isoIndex][8*(tin)+1] = pf1;
            pf[isoIndex][8*(tin)+2] = pf2;
            pf[isoIndex][8*(tin)+3] = pf3;
            pf[isoIndex][8*(tin)+4] = pf4;
            pf[isoIndex][8*(tin)+5] = pf5;
            pf[isoIndex][8*(tin)+6] = pf6;
            pf[isoIndex][8*(tin)+7] = pf7;
        }
        
        numberSpecies = isoIndex + 1;
    }
    
}    // End of function readNetwork (char *fileName)



// Function to compute the temperature-dependent ReacLib rates.  The factor
// Prefac contains the eta counting and sign parameter, and the corresponding
// mass density factors for each reaction. Thus computing fluxes in computeFluxCPU()
// only requires multiplying these rates by the corresponding abundance factors Y.

void computeRatesCPU(double T9)
{	
    double T93 = powf(T9, THIRD); 
    double t1 = 1/T9;
    double t2 = 1/T93;
    double t3 = T93;
    double t4 = T9;
    double t5 = T93*T93*T93*T93*T93;
    double t6 = logf(T9);
    
    for (int i=0; i<numberReactions; i++){
        Rate[i] = 
        Prefac[i]*expf(P0[i] + t1*P1[i] + t2*P2[i] + t3*P3[i] + t4*P4[i] + t5*P5[i] + t6*P6[i]);
    }
}    // End of function computeRatesCPU(double T9)



// Function to print out all the rates. The label can be used to distinguish cases
// if called more than once.

void  writeRates(char *label)
{
    printf("\n\nCOMPUTED RATES (%s):\n\n", label);
    for (int i=0; i<numberReactions; i++) 
    {
        printf("%d %s Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Q=%5.3f Prefac=%6.3e Reactants=%d\n",
               i,reacLabel[i], Rate[i], Y[Reactant1[i]], Y[Reactant2[i]], Y[Reactant3[i]], 
               Q[i], Prefac[i], NumReactingSpecies[i]);
    }
    
}    // End of function writeRates(char *label)


// Function to compute the current flux for each reaction

void computeFluxCPU()
{
    for(int j=0; j<numberReactions; j++)
    {
        int nr = NumReactingSpecies[j];
        
        switch(nr)
        {		
            case 1:
                Flux[j] = Rate[j] * Y[*(Reactant1+j)];	
                if(showFluxCalc == 1)
                {
                    printf("\nj=%d %s Rate=%7.3e Y1=%7.3e Flux=%7.3e",j,reacLabel[j],Rate[j],Y[*(Reactant1+j)],Flux[j]);
                }
                break;
                
            case 2:					
                Flux[j] = Rate[j] * Y[*(Reactant1+j)] * Y[*(Reactant2+j)]; 	
                if(showFluxCalc == 1)
                {
                    printf("\nj=%d %s Rate=%7.3e Y1=%7.3e Y2=%7.3e Flux=%7.3e",j,reacLabel[j],Rate[j],Y[*(Reactant1+j)],
                           Y[*(Reactant2+j)],Flux[j]);
                }
                break;
                
            case 3:					
                Flux[j] = Rate[j] * Y[*(Reactant1+j)] * Y[*(Reactant2+j)] 
                * Y[*(Reactant3+j)];
                if(showFluxCalc == 1)
                {
                    printf("\nj=%d %s Rate=%7.3e Y1=%7.3e Y2=%7.3e Y3=%7.3e Flux=%7.3e",j,reacLabel[j],Rate[j],Y[*(Reactant1+j)],
                           Y[*(Reactant2+j)],Y[*(Reactant3+j)],Flux[j]);
                }
                break;
        }
    }
    printf("\n");

}    // End of function computeFluxes()



// Function to print out the network isotopes, mass excesses, and the entries in the 
// partition function table for each isotope.

void writeNetwork()
{
    printf("\n\n%d ISOTOPES IN NETWORK:\n\n",numberSpecies);
    printf("Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    for (int i=0; i<numberSpecies; i++) 
    {
        printf("%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n",  i, isoLabel[i], (int)AA[i], Z[i], N[i], 
               Y[i], X[i], massExcess[i]);
    }
    
    printf("\nPARTITION FUNCTION TABLE:\n");
    printf("\n T9 = %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
    Tpf[0],Tpf[1],Tpf[2],Tpf[3],Tpf[4],Tpf[5],Tpf[6],Tpf[7],
    Tpf[8],Tpf[9],Tpf[10],Tpf[11],Tpf[12],Tpf[13],Tpf[14],Tpf[15],
    Tpf[16],Tpf[17],Tpf[18],Tpf[19],Tpf[20],Tpf[21],Tpf[22],Tpf[23]
    );
    for(int j=0; j<numberSpecies; j++){
        printf("\n%-5s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
        isoLabel[j],pf[j][0],pf[j][1],pf[j][2],pf[j][3],pf[j][4],pf[j][5],pf[j][6],pf[j][7],	
        pf[j][8],pf[j][9],pf[j][10],pf[j][11],pf[j][12],pf[j][13],pf[j][14],pf[j][15],
        pf[j][16],pf[j][17],pf[j][18],pf[j][19],pf[j][20],pf[j][21],pf[j][22],pf[j][23] );
    }
    printf("\n");
    
}   // End of function writeNetwork()



/*
 * Function to find the contributions to F+ and F- of each reaction for each isotope.  
 * This is executed only once at the beginning of the entire calculation to determine 
 * the structure of the network.
 */

void parseF()
{
    if(showParsing == 1)
        printf("\n--- Use parseF() to find F+ and F- flux components for each species ---");
    
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
        // change its population up (contributing to F+) or down (contributing
        // to F-).
        
        for(int j=0; j<numberReactions; j++)
        {
            int totalL = 0;
            int totalR = 0;
            
            // Loop over reactants for this reaction
            for(int k=0; k<NumReactingSpecies[j]; k++)
            {
                if(Z[i] == reactantZ[j][k] && N[i] == reactantN[j][k]) totalL ++;
            }
            
            // Loop over products for this reaction
            for(int k=0; k<NumProducts[j]; k++)
            {
                if(Z[i] == productZ[j][k] && N[i] == productN[j][k]) totalR ++;
            }
            
            total = totalL - totalR;
            
            if(total > 0)        // Contributes to F- for this isotope
            {
                numFminus ++;
                reacMask[i][j] = -total;
                tempInt2[incrementMinus + numFminus-1] = j;
                if(showParsing == 1)
                    printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F-", 
                           isoLabel[i], j, reacLabel[j], NumReactingSpecies[j], NumProducts[j], totalL, 
                           totalR, total);
            } 
            else if(total < 0)   // Contributes to F+ for this isotope
            {
                numFplus ++;
                reacMask[i][j] = -total;
                tempInt1[incrementPlus + numFplus-1] = j;
                if(showParsing == 1)
                    printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d totR=%d tot=%d F+", 
                           isoLabel[i], j, reacLabel[j], NumReactingSpecies[j], NumProducts[j], totalL, 
                           totalR, total);
            }
            else                 // Does not contribute to flux for this isotope
            {
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
            printf("\n%d %s numF+ = %d numF- = %d", i, isoLabel[i], numFplus, numFminus);
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


/* Use the information constructed in parseF (in particular reacMask[j][k]) to create the 
 reaction vectors for the network*/

void makeReactionVectors(){
    
    // Fill a 2D array RV[j][k] where j labels the reaction index and k labels
    // the components of the reaction vector for that reaction.
    
    for(int j=0; j<numberReactions; j++){
        for(int k=0; k<ISOTOPES; k++){
            RV[j][k] = reacMask[k][j];
            //printf("\n   j=%d k=%d RV=%d reacMask=%d", j, k, RV[j][k], reacMask[k][j]);
        }
    }
    
    // Write out the array containin components of the reaction vectors
    
    int uppity = minimumOf(25, numberSpecies);  // limit printout width to 25 species
    printf("\n\nREACTION VECTORS (%d Reaction vectors with %d species components):\n", 
        numberReactions, ISOTOPES);
    printf("\nReaction \\ Species           ");
    for(int k=0; k<uppity; k++){
        printf("%4s  ",isoLabel[k]);
    }
    printf("\n");
    for(int j=0; j<numberReactions; j++){
        printf("%4d %22s [ ",j,reacLabel[j]);
        for(int k=0; k<uppity-1; k++){
            printf("%2d    ",RV[j][k]);
        }
        printf("%2d ]\n",RV[j][uppity-1]);
    }
    
    
    // -----------------------------------------------------------------------
    // Now implement reaction vectors as GSL vectors so that we can use the
    // GSL and GSL_BLAS API to manipulate them.  The prototypes for doing this
    // may be found in the example code arrayPointers.c and matrixGSL.c
    // -----------------------------------------------------------------------
    
    rvPt = rv;             // Set pointer to beginning address of array rv
    
    // Allocate an array populated with GSL vectors
    
    printf("\nAllocating an array rv[] of GSL vectors\n");
    for(int i=0; i<SIZE; i++){
        
        //Prototype GSL reaction vector
        gsl_vector * v1 = gsl_vector_alloc (ISOTOPES); 
        
        // Set elements of rv[] pointed to by *rvPt equal to GSL vectors
        *(rvPt+i) = *v1;   
    }
    
    // Fill vector component entries created above with data contained in  
    // RV[i][j]
    
    printf("\nPopulate vector components of array rv[i]_j with RV[i][j]:");
    for (int i = 0; i < SIZE; i++)
    {
        printf("\n\nrv[%d]",i);
        for(int j=0; j<ISOTOPES; j++){
            gsl_vector_set (rvPt+i, j, RV[i][j]);
            // Retrieve the vector component just stored and print it
            int temp = gsl_vector_get(rvPt+i, j);
            printf("\ni=%d j=%d  RV[%d][%d] =%3d  rv[%d]_%d =%3d", 
                i, j, i, j, RV[i][j], i, j, temp);
        }
    }
    
//     // Access and print values of vector components populated above for each rv[] entry
//     
//     printf("\n\nAccess and print each component of array of GSL vectors rv[i]_j:");
//     for (int i = 0; i < SIZE; i++)
//     {
//         printf("\n");
//         for(int j=0; j<ISOTOPES; j++){
//             
//             // Define a pointer that will point to the GSL vector in array entry rv[i].
//             
//             gsl_vector *vector = rvPt+i;
//             
//             // Assign the jth component of the vector in rv[i] to a variable
//             
//             int component = gsl_vector_get (vector, j);
//             printf ("\nrv[%d]_%d = %d", i, j, component);
//         }
//     }
    
    // Display reaction vectors as component list
    
    printf("\n\nREACTION VECTOR COMPONENTS (%d reaction vectors with %d components)\n",
        SIZE, ISOTOPES);
    
    for (int i = 0; i < SIZE; i++)
    {
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
}


// ------------------------------------------------------------------------
// Method to use compareGSLvectors to sort all reactions in the network
// into reaction groups labeled by a series of integers 0, 1, ...  All 
// reactions in a reaction group have the same reaction vector up to a 
// sign. The array RGindex[] of dimension SIZE holds the integer labeling 
// reaction group for each reaction after this function is executed. 
// ------------------------------------------------------------------------

void sortReactionGroups(void){
    
    int showall = 0;    // 1 to show diagnostic printout; 0 to suppress
    
    // Cycle over all reaction vectors and compare them pairwise to 
    // assign to reaction groups. The pointer rvPt points to the array
    // rv[] of GSL reaction vectors.
    
    // The integer rindex labels the reaction group
    int rindex = -1;
    int ck = -1;
    
    // Initialize
    for(int i=0; i<SIZE; i++){
        RGindex[i] = -1;
        //printf("\ni=%d doneAlready=%d ",i, doneAlready[i]);
    }
    
    int scorekeeper = 0;
    for (int i=0; i<SIZE; i++){
        scorekeeper = 0;
        if(i==0) rindex ++;
        if(showall == 1) printf("\nRG=%d", rindex);
        for(int j=0; j<SIZE; j++){
            
            if(RGindex[i] < 0) RGindex[i] = rindex;
            ck = compareGSLvectors(rvPt+i, rvPt+j);
            
            //printf("\n%d Not done: %d", j, !doneAlready[j]);
            if(ck > 0 && RGindex[j]< 0) {
                RGindex[j] = rindex;
                scorekeeper ++;
                //doneAlready[j] = true;
            }
            if(showall==1)
            printf("\ni=%d j=%d RGindex[%d]=%d ck=%d rindex=%d scorekeeper=%d", 
                   i, j, j, RGindex[j], ck, rindex, scorekeeper);
            
            if(i>43)
            printf("\ni=%d j=%d numberMembers=%d rg=%d ck=%d RGindex[j]=%d",
                    i, j, scorekeeper+1, rindex, ck, RGindex[j]
            );
        }
        if(scorekeeper > 0) rindex++;
    }
    
    numberRG = rindex;   // Store total number of reaction groups
    
    // Diagnostic showing reaction group associated with each reaction
    if(showall == 1){
        for(int i=0; i<SIZE; i++){
            printf("\ni=%d  RGindex=%d", i, RGindex[i]);
        }
    }
    
    // Write out the components of the reaction groups
    
    printf("\n--- PARTIAL EQUILIBRIUM REACTION GROUPS ---");
    for(int i=0; i<numberRG; i++){
        printf("\n\nReaction group %d", i);
        for(int j=0; j<SIZE; j++){
            if(RGindex[j] == i){
                printf("\n  %d %s",  j, reacLabel[j]);
            }
        }
    }
    printf("\n");
    
}      // End sortReactionGroups()


// ------------------------------------------------------------------------
// Method to compare two GSL vectors of same length.  Returns 0 
// if they are not equivalent, 1 if they are the same, 2 if one vector is 
// the negative of the other. The arguments are pointers to the two
// GSL vectors.
// ------------------------------------------------------------------------

int compareGSLvectors(gsl_vector* rv1, gsl_vector* rv2){
    
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
}


// -------------------------------------------------------------------------
// Utility function to return the network vector index given Z and N.
// -------------------------------------------------------------------------

int returnNetIndexZN(int z, int n) {
    for (int i = 0; i < numberSpecies; i++) {
        if (Z[i] == z && N[i] == n)
            return i;
    }
    return -1;
}


// -----------------------------------------------------------------------
// Utility function to return the network vector index given the symbol.
// -----------------------------------------------------------------------

int returnNetIndexSymbol(char* symbol) {
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
// Utility function to return true if given (Z,N) in the network, false 
// otherwise.
// ----------------------------------------------------------------------

bool isInNet(int Z, int N) {
    if (returnNetIndexZN(Z, N) < 0) {
        return false;
    } else {
        return true;
    }
}

// ----------------------------------------------------------------------
// Utility function to return minimum of two integers 
// ----------------------------------------------------------------------

int minimumOf(int x, int y) 
{ 
    return (x < y) ? x : y;
}

// ----------------------------------------------------------------------
// Utility function to return maximum of two integers 
// ----------------------------------------------------------------------

int maximumOf(int x, int y) 
{ 
    return (x > y) ? x : y; 
}


// Function to test the CPU timer by executing a long, pointless loop.
void testTimerCPU()
{
    double a, b;
    
    START_CPU;
    for (long count = 1l; count < 500000l; count++) {
        a = sqrt(count);
        b = 1.0/logf(a);
        a = logf(b)/sqrt(a);
    }
    STOP_CPU;
    PRINT_CPU;
    
}   // End of function testTimerCPU()
