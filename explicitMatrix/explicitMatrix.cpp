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
//        filenames rateLibraryFile[] and networkFile[] to appropriate values.

// SIZE defines the number of reactions to be calculated.  Need min of SIZE=4395 for 365-element network, 
// 1603 for 150-isotope network, 2234 for the 194-isotope network, 3176 for the 268-isotope network,
// 1566 for the nova134 network, 48 for the alpha network, 8 for 3-alpha 4He-12C-16O network, and 28 for 
// the 7-isotope pp network. ISOTOPES defines the number of isotopes in each network (for example, ISOTOPES=7
// for the 7-isotope pp-chain network.  These sizes are hardwired for now but eventually we may want to read
// them in and assign them dynamically.

#define ISOTOPES 3                    // Max isotopes in network (e.g. 16 for alpha network)
#define SIZE 8                        // Max number of reactions (e.g. 48 for alpha network)

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
#define PRINT_CPU (printf("\nTimer: %g ms used by CPU\n", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
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
int displayInput = 1;
int showParsing = 1;
int showFparsing = 1;
int showFluxCalc = 1;

// Function Signatures:
void devcheck(int);
void readLibraryParams(char *);
void readNetwork(char *);
void writeNetwork(void);
void writeRates(char *);
void writeAbundances(void);
void setRG(int, int, int);

double Rate[SIZE];       // Computed rate for each reaction from Reactions::computeRate()
double Flux[SIZE];       // Computed flux for each reaction from Reactions::computeFlux()

// ------ Species data in the following arrays also contained in class Species

int Z[ISOTOPES];                 // Array holding Z values for isotopes
int N[ISOTOPES];                 // Array holding N values for isotopes
int AA[ISOTOPES];             // Array holding A values for isotopes
double Y[ISOTOPES];              // Array holding abundances Y for isotopes
double X[ISOTOPES];              // Array holding mass fractions X for isotopes
double massExcess[ISOTOPES];     // Array holding mass excesses for isotopes
char isoLabel[ISOTOPES][5];      // Isotope labels (max length 5 characters; e.g. 238pu)

// -----------


// ------ Reaction data in the following arrays also contained in class Reaction

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

int numberRG;          // Number of partial equilibrium reaction groups

// Define array to hold the reaction group index for each reaction. There are n reaction
// groups in the network and each reaction belongs to one reaction group.  RGindex[m] is
// the index (0, 1, ... n) of the reaction group to which reaction m belongs. 
// This array is populated by the function ReactionGroup::sortReactionGroups()

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


/* Class Utilities to hold utility useful utility functions.  Functions are
 * declared static so that they can be invoked without having to instantiate
 * objects of type Utilities.  For example, Utilities::returnNetIndexZN (Z, N).
 * We will generally have other classes inherit from Utilities so that they can
 * access its static methods.  This class definition must precede definitions of 
 * classes that inherit from it.
 */


class Utilities{
    
private:
    
public:
    
    // -------------------------------------------------------------------------
    // Static function to return the network vector index given Z and N
    // for the isotope.  Returns -1 if no match.
    // -------------------------------------------------------------------------
    
    static int returnNetIndexZN(int z, int n) {
        for (int i = 0; i < numberSpecies; i++) {
            if (Z[i] == z && N[i] == n)
                return i;
        }
        return -1;
    }
    
    
    // -----------------------------------------------------------------------
    // Static function to return the network vector index given the symbol
    // for the isotope.
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
    // Static function to return true if given (Z,N) in the network, false 
    // otherwise.
    // ----------------------------------------------------------------------
    
    static bool isInNet(int Z, int N) {
        if (returnNetIndexZN(Z, N) < 0) {
            return false;
        } else {
            return true;
        }
    }
    
    // ----------------------------------------------------------------------
    // Static function to return minimum of two numbers.  Overload to 
    // accept either integer or double arguments.
    // ----------------------------------------------------------------------
    
    static int minimumOf(int x, int y) {   
        return (x < y) ? x : y;
    }
    
    static double minimumOf(double x, double y) { 
        return (x < y) ? x : y;
    }
    
    // ----------------------------------------------------------------------
    // Static function to return maximum of two numbers.  Overload to 
    // accept either integer or double arguments.
    // ----------------------------------------------------------------------
    
    static int maximumOf(int x, int y) { 
        return (x > y) ? x : y; 
    }
    
    static double maximumOf(double x, double y){ 
        return (x > y) ? x : y; 
    }
    
    // ----------------------------------------------------------------------
    // Static function to convert string to char so that it will print
    // in printf. Assumes max of 50 characters.  Typically a string can
    // be printed with cout but a string given to printf typically displays
    // garbage because of type issues in the C function printf. This
    // function converts a string to a corresponding character array, which
    // printf or cout can print.
    // ----------------------------------------------------------------------
    
    static char* stringToChar(string s){
        static char cs[50];
        strcpy(cs, &s[0]);   // alternatively strcpy(cs, s.c_str());
        return cs;
    }
    
    // ----------------------------------------------------------------------
    // Static unction to test CPU timer by executing long, pointless loop.
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
 * for each isotope in the network. Inherits from class Utilities.
 */

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
        
        
};      // End of class Species



/* Class Reaction to describe all reactions in the network.  Create new instance
 * for each reaction in the network. Inherits from Utilities
 */

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
            //int i; 
            for (int i = 0; i < sizeof(p); i++) { 
                p[i] = s[i]; 
                reacLabel[reacIndex][i] = p[i];
            }
        }
        
        void setisEC(int i){ isEC = i; }
        
        void setisReverse(int i){ isReverse = i; }
        
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
        
        double getQ(){ return Q; }
        
        double getprefac(){ return prefac; }
        
        double getp(int k){ return p[k]; }
        
        int getreactantZ(int k){
            if(k > numberReactants){
                printf("\n\nERROR: k=%d larger than number reactants %d", k, numberReactants);
                return -1;
            } else {
                return reactantZ[k];
            }
        }
        
        int getreactantN(int k){
            if(k > numberReactants){
                printf("\n\nERROR: k=%d larger than number reactants %d", k, numberReactants);
                return -1;
            } else {
                return reactantN[k];
            }
        }
        
        int getproductZ(int k){
            if(k > numberProducts){
                printf("\n\nERROR: k=%d larger than number products %d", k, numberProducts);
                return -1;
            } else {
                return productZ[k];
            }
        }
        
        int getproductN(int k){
            if(k > numberProducts){
                printf("\n\nERROR: k=%d larger than number products %d", k, numberProducts);
                return -1;
            } else {
                return productN[k];
            }
        }
        
        int getreactantIndex(int k){
            if(k > numberReactants){
                printf("\n\nERROR: k=%d larger than number reactants %d", k, numberReactants);
                return -1;
            } else {
                return reactantIndex[k];
            }
        }
        
        int getproductIndex(int k){
            if(k > numberProducts){
                printf("\n\nERROR: k=%d larger than number products %d", k, numberProducts);
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
            
            printf("\n-- Initial values of F+\n\n");
            for(int i=0; i<totalFplus; i++){
                int indy = MapFplus[i];
                Fplus[i] = FplusFac[i]*Flux[indy];
                printf("i=%d FplusFac=%3.1f Flux=%7.3e Fplus=%7.3e\n", 
                       i, FplusFac[i], Flux[indy], Fplus[i]);
            }
            
            printf("\n\n-- Initial values of F-\n\n");
            
            for(int i=0; i<totalFminus; i++){
                int indy = MapFminus[i];
                Fminus[i] = FminusFac[i]*Flux[indy];
                printf("i=%d FminusFac=%3.1f Flux=%7.3e Fminus=%7.3e\n", 
                       i, FminusFac[i], Flux[indy], Fminus[i]);
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
            
            printf("\nReaction::computeRate  %d %19s densfac=%6.3e rate= %6.3e Rrate=%6.3e", 
                   getreacIndex(), getreacChar(), getdensfac(), getrate(), getRrate());
            
        }
        
        
        // Reaction::computeFlux() to compute the current flux for reaction corresponding 
        // to this Reaction object.
        
        void computeFlux(){
            
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
                        printf("\n%d %18s reactants=%d iso0=%d iso1=%d Rrate=%7.3e Y1=%7.3e Y2=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], reactantIndex[1], Rrate,
                            Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], flux);
                    }
                    break;
                    
                case 3:	   // 3-body reactions
                    
                    flux = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ] * Y[ reactantIndex[2] ];
                    Flux[getreacIndex()] = flux;         // Put in flux array in main
                    if(showFluxCalc == 1){
                        printf("\n%d %18s reactants=%d iso0=%d iso1=%d iso2=%d Rrate=%7.3e Y1=%7.3e Y2=%7.3e Y3=%7.3e Flux=%7.3e",
                            reacIndex, getreacChar(), numberReactants, reactantIndex[0], reactantIndex[1], reactantIndex[2], 
                            Rrate, Y[ reactantIndex[0] ], Y[ reactantIndex[1] ], 
                            Y[ reactantIndex[2] ], flux);
                    }
                    break;
            }
            
        }    // End of function computeFlux()
        
        
};  // End class Reaction



// Class ReactionVector with  methods that create reaction vectors. Inherits from 
// class Utilities.


class ReactionVector:  public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
private:
    
    // 1 to print diagnostic, 0 to suppress
    
    static const int printExtra = 0;
    static const int showall = 1;    
    
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
        
        printf("\nAllocating an array rv[] of GSL vectors\n");
        for(int i=0; i<SIZE; i++){
            
            //Prototype GSL reaction vector
            gsl_vector * v1 = gsl_vector_alloc (ISOTOPES); 
            
            // Set elements of rv[] pointed to by *rvPt equal to GSL vectors
            *(rvPt+i) = *v1;   
        }
        
        // Fill vector component entries created above with data contained in  
        // reacMask[j][i] (notice reversed indices)
        
        if(printExtra == 1){
            printf("\nPopulate vector components of array rv[i]_j with reacMask[j][i]:");
        }
        
        for (int i = 0; i < SIZE; i++) {
            if(printExtra == 1) printf("\n\nrv[%d]",i);
            for(int j=0; j<ISOTOPES; j++){
                gsl_vector_set (rvPt+i, j, reacMask[j][i]);
                // Retrieve the vector component just stored and print it
                int temp = gsl_vector_get(rvPt+i, j);
                if(printExtra == 1){
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
        
        printf("\n\n\n--- SORTING REACTION GROUPS ---");
        
        int scorekeeper = 0;
        for (int i=0; i<SIZE; i++){
            scorekeeper = 0;
            if(i==0) rindex ++;
            if(showall == 1) printf("\n\nRG=%d", rindex);
            for(int j=0; j<SIZE; j++){
                
                if(RGindex[i] < 0) RGindex[i] = rindex;
                ck = compareGSLvectors(rvPt+i, rvPt+j);
                
                if(ck > 0 && RGindex[j]< 0) {
                    RGindex[j] = rindex;
                    scorekeeper ++;
                }
                if(showall==1){
                    printf("\ni=%d j=%d RGindex[%d]=%d ck=%d rindex=%d scorekeeper=%d", 
                        i, j, j, RGindex[j], ck, rindex, scorekeeper);
                }
            }
            
            if(scorekeeper > 0) rindex++;
        }
        
        numberRG = rindex;   // Store total number of reaction groups
        
        // Diagnostic showing reaction group associated with each reaction
        
        if(showall == 1){
            printf("\n\n-- SUMMARY OF REACTION GROUPS:\n");
            for(int i=0; i<SIZE; i++){
                printf("\nreaction=%d  %18s RGindex=%d RGmemberIndex=%d", 
                       i, reacLabel[i], RGindex[i], RGMemberIndex[i]);
            }
        }
        
        // Write out the components of the reaction groups
        
        printf("\n\n\n--- PARTIAL EQUILIBRIUM REACTION GROUPS ---");
        for(int i=0; i<numberRG; i++){
            printf("\n\n-- Reaction Group %d", i);
            int rgindex = -1;
            for(int j=0; j<SIZE; j++){
                if(RGindex[j] == i){
                    rgindex ++; 
                    setRG(j, RGclass[j], RGindex[j]);
                    printf("\n%s reacIndex=%d RGindex=%d RG=%d RGreacIndex=%d RG: %s", 
                           reacLabel[j], j, rgindex, RGclass[j], RGMemberIndex[j],
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




// Main CPU routine

int main() { 
    
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
    
    // Read in network file and associated partition functions.  This is required only
    // once at the beginning of the entire calculation.  
    
    char *networkFilePtr = networkFile;
    readNetwork(networkFilePtr);
    writeNetwork();
    
    // Read in rate library data from a file. This is required only once, at the
    // beginning of the entire calculation.
    
    char *rateLibraryFilePtr = rateLibraryFile;
    readLibraryParams(rateLibraryFilePtr);
    
    // Print out some quantitites from the Reaction object reaction[]
    
    // First, print using the string class object reacString with cout
    
    printf("\n\nQuantities stored in Reaction class objects reaction[]\n\n");
    printf("Print with cout without converting strings to chars:\n");
    for(int i=0; i<SIZE; i++){
        cout << reaction[i].getreacIndex() << " "
             << reaction[i].getreacString() 
             << " reacClass=" << reaction[i].getreacClass()
             << " reactants=" << reaction[i].getnumberReactants()
             << " products=" << reaction[i].getnumberProducts()
             << " isEC=" << reaction[i].getisEC()
             << " isReverse=" << reaction[i].getisReverse()
             << " Q=" << reaction[i].getQ()
             << " prefac=" << reaction[i].getprefac()
             << endl ;
    }
    
    // Alternatively, convert the string class object reacString to a char and
    // print with printf.
    
    printf("\n");
    printf("\nAlternatively, formatted print with printf after converting strings to chars:\n");
    string ts = reaction[0].getreacString();
    char tsc[ts.size() + 1];     // Allocate temporary character array
    
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
    
    printf("\n\nReacLib Parameters read in for the %d network reactions:\n", SIZE);
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
            printf(" reactantIndex[%d]=%d", j, reaction[i].getproductIndex(j));
        }
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
    printf("\n\nNETWORK SPECIES VECTOR (%d components):\n\nIndex  Species    Z     N",numberSpecies);
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
    
    // Use methods of Reaction class to compute reaction rates. We have instantiated
    // a set of Reaction objects in the array reaction[i], one entry for each
    // reaction in the network. Loop over this array and call the computeRate()
    // method of Reaction on each object.
    
    printf("\nCOMPUTED RATES using Reaction::computeRate(T9, rho)\n");
    
    START_CPU     // Start a timer for rate calculation
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeConstantFacs(T9, rho);
        reaction[i].computeRate(T9, rho);
    }
    
    STOP_CPU;     // Stop the timer
    printf("\n");
    PRINT_CPU;    // Print timing information for rate calculation
    printf("\n");
    
    // Use methods of the Reaction class to compute fluxes.  We have instantiated
    // a set of Reaction objects in the array reaction[i], one entry for each
    // reaction in the network. Loop over this array and call the computeFlux()
    // method of Reaction on each object.
    
    printf("\n-- TOTAL FLUXES (T9=%4.2f, rho=%7.3e) from Reaction::computeFlux()\n", 
        T9, rho);
    
    START_CPU     // Start a timer for rate calculation
    
    for(int i=0; i<SIZE; i++){
        reaction[i].computeFlux();
    }
    
    STOP_CPU;        // Stop the timer
    printf("\n");
    PRINT_CPU;       // Print timing information for flux calculation
    printf("\n");
    
    // Call the static function Reaction::populateFplusFminus() to populate F+ and F-
    // for each isotope set up in setupFplusFminus() from master flux array.
    
    Reaction::populateFplusFminus();
  
    
    // Perform some tests of various functions
    
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
    // using pointers instead.
    
    // Set Species pointer to address of Species object isotope[1]
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
    //int check = reactionVector.compareGSLvectors(rvPt+indy1, rvPt+indy2);
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
    
    // Test of CPU timer by executing a long, pointless loop
    Utilities::testTimerCPU();

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
    int i0, i1, i2, i3, i4, i5, i6;
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
                sscanf (line, "%s %d %d %d %d %d %d %d %lf %lf", rlabel, &i0, &i1, &i2, &i3, &i4, &i5, &i6, 
                        &sf, &q);
                
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
    
    printf("\n\nnumberReactions=%d",numberReactions);
    
    fclose(fr);           // Close the file
    
}    // End of function readLibraryParams (char *fileName)



// Function to print out the network isotopes, mass excesses, and the entries in the 
// partition function table for each isotope.

void writeNetwork()
{
    printf("\n\n%d ISOTOPES IN NETWORK:\n\n",numberSpecies);
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
