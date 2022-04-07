#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// NOTE:  For different networks you must change the values of SIZE, ISOTOPES, rateLibraryFile[],
//        and networkFile[] to appropriate values.

// SIZE defines the number of reactions to be calculated.  Need min of SIZE=4395 for 365-element network, 
// 1603 for 150-isotope network, 2234 for the 194-isotope network, 3176 for the 268-isotope network,
// 1566 for the nova134 network, 48 for the alpha network, 8 for 3-alpha 4He-12C-16O network, and 28 for 
// the 7-isotope pp network. ISOTOPES defines the number of isotopes in each network (for example, SIZE=7
// for the 7-isotope pp-chain network.  These sizes are hardwired for now but eventually we want to read
// them in and assign them dynamically.

#define SIZE 28                       // Max number of reactions (48 for alpha network)
#define ISOTOPES 7                    // Max isotopes in network (16 for alpha network) 16

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

char rateLibraryFile[] = "rateLibrary_pp.data"; //"rateLibrary_3alpha.data"; //"rateLibrary_alpha.data";

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format for 
// this file. Standard test cases: CUDAnet_alphasolar.inp, CUDAnet_150solar.inp,
// CUDAnet_365solar.inp, CUDAnet_nova134.inp, CUDAnet_3alpha.inp, CUDAnet_pp.inp.

char networkFile[] = "CUDAnet_pp.inp"; //"CUDAnet_3alpha.inp"; //"CUDAnet_alpha.inp";

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
void computeRatesCPU(float);
void computeFluxCPU(void);
void writeRates(char *);
void writeAbundances(void);
void parseF(void);

// Variables to hold data that will be read in

float P0[SIZE];       // Array holding library rate parameters p0
float P1[SIZE];       // Array holding library rate parameters p1
float P2[SIZE];       // Array holding library rate parameters p2
float P3[SIZE];       // Array holding library rate parameters p3
float P4[SIZE];       // Array holding library rate parameters p4
float P5[SIZE];       // Array holding library rate parameters p5
float P6[SIZE];       // Array holding library rate parameters p6
float Q[SIZE];        // Array holding library entry for reaction Q-value
float Rate[SIZE];     // Array holding the computed rate for each reaction
float Flux[SIZE];     // Array holding the computed fluxes for each reaction

int reaclibClass[SIZE];         // Reaclib class index for reaction (1-8) 
int RGclass[SIZE];              // Reaction Group class (PE) for reaction (1-5)
int RGmemberIndex[SIZE];        // Member index within reaction group
int NumReactingSpecies[SIZE];         // Number of reactant isotopes for the reaction
int NumProducts[SIZE];          // Number of product isotopes for the reaction
int isEC[SIZE];                 // Whether reaction is electron capture (0=false; 1=true)
int isReverseR[SIZE];           // Whether reaction is inverse (matters for partition func)
float Prefac[SIZE];             // The statistical prefactor for each reaction

int reactantZ[SIZE][4];          // Holds Z for each reactant isotope
int reactantN[SIZE][4];          // Holds N for each reactant isotope
int productZ[SIZE][4];           // Holds Z for each product isotope
int productN[SIZE][4];           // Holds N for each product isotope
int ReactantIndex[SIZE][4];      // Index of isotope vector for each reactant isotope
int ProductIndex[SIZE][4];       // Index of isotope vector for each product isotope
int Reactant1[SIZE];             // Isotope index of first reactant
int Reactant2[SIZE];             // Isotope index of second reactant (if 2-body or 3-body)
int Reactant3[SIZE];		     // Isotope index of third reactant (if 3-body)

int Z[ISOTOPES];                // Array holding Z values for isotopes
int N[ISOTOPES];                // Array holding N values for isotopes
float AA[ISOTOPES];             // Array holding A values for isotopes
float Y[ISOTOPES];              // Array holding abundances Y for isotopes
float X[ISOTOPES];              // Array holding mass fractions X for isotopes
float massExcess[ISOTOPES];     // Array holding mass excesses for isotopes
char isoLabel[ISOTOPES][5];     // Isotope labels

int numberSpecies;              // Total number of species in network
int numberReactions;            // Total number of reactions in the network

char reacLabel[SIZE][LABELSIZE];   // Character array holding reaction labels

// Temperatures in units of 10^9 K for partition function table (see pf[][]). These are
// copied from the corresponding variable Tpf[] in the Java code.

const float Tpf[] = { 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f,
    1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f };
    
// Array holding partition function values for the 24 temperatures given in Tpf[]
// for each isotope

float pf[ISOTOPES][PF]; 

// Dens[] will hold the density factor for fluxes. Dens[0] = 1.0 (1-body reactions), 
// Dens[1]=rho (2-body reactions), Dens[2]=rho*rho (3-body reactions)

float Dens[3];  

// Array with entries +1 if a reaction increases the population of the isotope (contributes to 
// F+), -1 if it decreases it (contributes to F-) and 0 if the reaction does not change the population
// of the isotope. This array is populated in the function parseF().  It is characteristic of the
// structure of the network and thus has to be calculated only once for a given network.

int reacMask[ISOTOPES][SIZE];

// Total number of F+ and F- terms in the network
int totalFplus = 0;
int totalFminus = 0;

// Arrays to hold non-zero fluxes in the network. Corresponding memory will be allocated dynamically 
// below with malloc

float* Fplus;        // Dynamically-allocated 1D array for non-zero F+ (Dim totalFplus)
float* Fminus;       // Dynamically-allocated 1D array for non-zero F- (Dim totalFminus)

// Arrays to hold number species factors for F+ and F- arrays. For example, for 12C+12C ->
// 4He + 20Ne the species factor is 1 for 4He and 20Ne diff. equation terms but 2 for
// 12C diff equation term (since the reaction involves one 4He and one 20Ne, but two 12C.

float* FplusFac;     // Dynamically allocated 1D array for number species factor in F+ terms
float* FminusFac;    // Dynamically allocated 1D array for number species factor in F- terms

float* FplusSum;     // Sum of F+ for each isotope
float* FminusSum;    // Sum of F- for each isotope

int* FplusMax;       // Upper index for each isotope in the Fplus array
int* FplusMin;       // Lower index for each isotope in the Fplus array
int* FminusMax;      // Upper index for each isotope in the Fminus array
int* FminusMin;      // Lower index for each isotope in the Fminus array

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


// Main CPU routine

int main()
{  
    
    // Set the temperature in units of 10^9 K and density in units of g/cm^3. In a
    // realistic calculation the temperature and density will be passed from the hydro 
    // code in an operator-split coupling of this network to hydro. Here we hardwire
    // them for testing purposes.  These will be used to calculate the reaction
    // rates in the network. Since we are assuming operator splitting, the temperature 
    // and density are assumed constant for each network integration.
    
    float T9 = 6.0f;
    float rho = 1.0e8;
    
    // Set the range of time integration and the initial timestep.  In an operator-split
    // coupling tmax will come from the hydro and dt_init will likely be the last timestep
    // of the previous network integration (for the preceding hydro timestep). Here we
    // hardwire them for testing purposes.
    
    float tmax = 1e-11;
    float dt_init = 1e-17; 
    
    // Read in network file and associated partition functions.  This is required only
    // once at the very beginning of the calculation.  
    
    char *networkFilePtr = networkFile;
    readNetwork(networkFilePtr);
    writeNetwork();
    
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
        
//         // Test the CPU timer by executing a long, pointless loop
//         testTimerCPU();
         
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
    
    parseF();
    
    // Allocate 1D arrays to hold non-zero F+ and F- for all reactions for all isotopes,
    // the arrays holding the species factors FplusFac and FminusFac, 
    // and also arrays to hold their sums for each isotope. Note that parseF() must
    // be run first because it determines totalFplus and totalFminus.
    
    Fplus = (float*) malloc(sizeof(float) * totalFplus);
    Fminus = (float*) malloc(sizeof(float) * totalFminus);
    FplusFac = (float*) malloc(sizeof(float) *totalFplus);
    FminusFac = (float*) malloc(sizeof(float) * totalFminus);
    FplusSum = (float*) malloc(sizeof(float) * numberSpecies);
    FminusSum = (float*) malloc(sizeof(float) * numberSpecies);
    
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
                FplusFac[tempCountPlus] = (float)reacMask[i][j];
// 				printf("\n F+  tempCountPlus=%d i=%d j=%d FplusFac=%3.1f", tempCountPlus, 
// 					   i, j, FplusFac[tempCountPlus]);
                tempCountPlus ++;
            }
            else if(reacMask[i][j] < 0)
            {
                FminusFac[tempCountMinus] = -(float) reacMask[i][j];
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
    
    printf("\n--- Values of F+ ---\n\n");
    for(int i=0; i<totalFplus; i++)
    {
        int indy = MapFplus[i];
        Fplus[i] = FplusFac[i]*Flux[indy];
        printf("i=%d FplusFac=%3.1f Flux=%7.3e Fplus=%7.3e\n", i,FplusFac[i],Flux[indy],Fplus[i]);
    }
    
    printf("\n--- Values of F- ---\n\n");
    
    for(int i=0; i<totalFminus; i++)
    {
        int indy = MapFminus[i];
        Fminus[i] = FminusFac[i]*Flux[indy];
        printf("i=%d FminusFac=%3.1f Flux=%7.3e Fminus=%7.3e\n", i,FminusFac[i],Flux[indy],Fminus[i]);
    }
    
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
    float p0, p1, p2, p3, p4, p5, p6, q, sf;
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
     *     float float float float float float float string
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
                sscanf (line, "%s %d %d %d %d %d %d %d %f %f", reaction, &i0, &i1, &i2, &i3, &i4, &i5, &i6, 
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
                sscanf (line, "%f %f %f %f %f %f %f", &p0, &p1, &p2, &p3, &p4, &p5, &p6);
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
 *	string int int int float float
 *	float float float float float float float float
 *	float float float float float float float float
 *	float float float float float float float float
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
    float y, mass;
    float pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7;
    
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
    
    if(displayInput==1) printf("\n\n--- Read in network and partition function ---\n");
    
    // Read in lines until NULL encountered. Lines can contain up to 60 characters
    
    while(fgets(line, 60, fr) != NULL)
    {
        isoSubIndex ++;
        if(isoSubIndex == 4){
            isoSubIndex = 0;
            isoIndex ++;
            // Scan and parse a title line
            sscanf (line, "%s %d %d %d %f %f", isoSymbol, &a, &z, &n, &y, &mass);
            if(displayInput == 1)
            {
                printf("\n%s %d %d %d %f %f\n", isoSymbol, a, z, n, y, mass);
            }
            // Store variables in arrays
            Z[isoIndex] = z;
            N[isoIndex] = n;
            AA[isoIndex] = (float)a;
            Y[isoIndex] = y;
            X[isoIndex] = AA[isoIndex]*Y[isoIndex];
            massExcess[isoIndex] = mass;
            for(int j=0; j<5; j++){
                isoLabel[isoIndex][j] = isoSymbol[j];
            }
        } else {
            // Scan and parse a partition function line. 
            sscanf (line, "%f %f %f %f %f %f %f %f", &pf0, &pf1, &pf2, &pf3, &pf4, &pf5, &pf6, &pf7);
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

void computeRatesCPU(float T9)
{	
    float T93 = powf(T9, THIRD); 
    float t1 = 1/T9;
    float t2 = 1/T93;
    float t3 = T93;
    float t4 = T9;
    float t5 = T93*T93*T93*T93*T93;
    float t6 = logf(T9);
    
    for (int i=0; i<numberReactions; i++){
        Rate[i] = 
        Prefac[i]*expf(P0[i] + t1*P1[i] + t2*P2[i] + t3*P3[i] + t4*P4[i] + t5*P5[i] + t6*P6[i]);
    }
}    // End of function computeRatesCPU(float T9)



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
    
    // Display some cases
    
    
    if(numberSpecies != 16 && numberSpecies > 25)
    {
        // Comment out this part of the if-block for alpha network to prevent warnings about index being
        // out of bounds. (Doesn't matter in calculation since this block is not reached if it is an alpha
        // network with 16 species, but the compile generates a long string of warnings.)  Uncomment
        // to show up to 26 species for larger networks.
        
        // 		printf("\n\nIndex               Reaction%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s\
        // %5s%5s%5s%5s%5s%5s",
        // 			isoLabel[0], isoLabel[1], isoLabel[2], isoLabel[3], isoLabel[4], isoLabel[5], isoLabel[6],
        // 			isoLabel[7], isoLabel[8], isoLabel[9], isoLabel[10], isoLabel[11], isoLabel[12], isoLabel[13],
        // 			isoLabel[14], isoLabel[15], isoLabel[16], isoLabel[17], isoLabel[18], isoLabel[19],
        // 			isoLabel[20], isoLabel[21], isoLabel[22], isoLabel[23], isoLabel[24], isoLabel[25]
        // 		);
        // 		for(int j=0; j<numberReactions; j++)
        // 		{
        // 			
        // 			printf(
        // 				"\n %4d %22s %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\
        //  %4d %4d %4d %4d %4d %4d",
        // 				j, reacLabel[j], reacMask[0][j], reacMask[1][j], reacMask[2][j],
        // 				reacMask[3][j], reacMask[4][j], reacMask[5][j], reacMask[6][j], 
        // 				reacMask[7][j], reacMask[8][j], reacMask[9][j], reacMask[10][j],
        // 				reacMask[11][j], reacMask[12][j], reacMask[13][j], reacMask[14][j],
        // 				reacMask[15][j], reacMask[16][j], reacMask[17][j], reacMask[18][j],
        // 				reacMask[19][j], reacMask[20][j], reacMask[21][j], reacMask[22][j],
        // 				reacMask[23][j], reacMask[24][j], reacMask[25][j]
        // 			);
        // 		}
        
    } 
    else if(numberSpecies > 15)  // For alpha networks
    {
        printf("\n\nFLUX-ISOTOPE COMPONENT ARRAY (-n --> F-; +n --> F+ for given isotope):");
        printf("\n\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
        printf("\n\nIndex               Reaction%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s",
               isoLabel[0], isoLabel[1], isoLabel[2], isoLabel[3], isoLabel[4], isoLabel[5], isoLabel[6],
               isoLabel[7], isoLabel[8], isoLabel[9], isoLabel[10], isoLabel[11], isoLabel[12], isoLabel[13],
               isoLabel[14], isoLabel[15]	
        );
        for(int j=0; j<numberReactions; j++)
        {
            
            printf(
                "\n %4d %22s %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d",
                j, reacLabel[j], reacMask[0][j], reacMask[1][j], reacMask[2][j],
                reacMask[3][j], reacMask[4][j], reacMask[5][j], reacMask[6][j], 
                reacMask[7][j], reacMask[8][j], reacMask[9][j], reacMask[10][j],
                reacMask[11][j], reacMask[12][j], reacMask[13][j], reacMask[14][j],
                reacMask[15][j]
            );
        }
        printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.\n\n", 
               totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
    }
    else if(numberSpecies==3)  // For 4Hi-12C-16O network
    {
        printf("\n\nFLUX-ISOTOPE COMPONENT ARRAY (-n --> F-; +n --> F+ for given isotope):");
        printf("\n\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
        printf("\n\nIndex               Reaction%5s%5s%5s",
               isoLabel[0], isoLabel[1], isoLabel[2]	
        );
        for(int j=0; j<numberReactions; j++)
        {
            
            printf(
                "\n %4d %22s %4d %4d %4d",
                j, reacLabel[j], reacMask[0][j], reacMask[1][j], reacMask[2][j]
            );
        }
        printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.\n", 
               totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
    }
    else if(numberSpecies==7)  // For pp chains
    {
        printf("\n\nFLUX-ISOTOPE COMPONENT ARRAY (-n --> F-; +n --> F+ for given isotope):");
        printf("\n\nnumberSpecies=%d numberReactions=%d",numberSpecies,numberReactions);
        printf("\n\nIndex               Reaction%5s%5s%5s%5s%5s%5s%5s",
               isoLabel[0], isoLabel[1], isoLabel[2], isoLabel[3], isoLabel[4], isoLabel[5], isoLabel[6]	
        );
        for(int j=0; j<numberReactions; j++)
        {
            
            printf(
                "\n %4d %22s %4d %4d %4d %4d %4d %4d %4d",
                j, reacLabel[j], reacMask[0][j], reacMask[1][j], reacMask[2][j],
                reacMask[3][j], reacMask[4][j], reacMask[5][j], reacMask[6][j]
            );
        }
        printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.\n", 
               totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
    }
    
}   // End of function parseF()



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
