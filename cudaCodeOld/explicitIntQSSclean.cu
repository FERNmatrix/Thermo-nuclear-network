
/* Program to calculate explicit QSS or asymptotic integration of a thermonuclear network
 * over a hydrodynamical timestep (temperature and density assumed constant over hydro timestep)
 * using a GPU. The goal is to do the entire network calculation on the GPU by copying the hydro
 * info (temperature, density) and abundances at the end of the hydro timestep from the CPU to
 * the GPU, running the network entirely on the GPU for the duration of the hydro timestep, and then
 * returning the abundances and energy release to the CPU for transmitting to the hydro for the next 
 * timestep, leaving the network intact on the GPU ready for the integration over the next hydro
 * timestep.

Mike Guidry
guidry@utk.edu
July 5, 2012
 
*/

#define BLOCKSIZE 512
#define DIAGNOSE_SIZE 300

// Kernel and device function definitions.  These will run on the device (the GPU).
// There is a max number of arguments and we are close to it. Adding 2-3 more puts us over.

__global__ void integrateNetwork (
	float* P0, 
	float* P1, 
	float* P2, 
	float* P3, 
	float* P4, 
	float* P5, 
	float* P6, 
	float* Prefac,  
	float* Q, 
	float* Rate, 
    float* Flux,
	float* Fplus,
	float* Fminus,
	float* FplusFac,
	float* FminusFac,
	float* FplusSum,
	float* FminusSum, 
 	int* FplusMax,
 	int* FminusMax,
	int* MapFplus,
	int* MapFminus,
	float* Y, 
	float* Diagnose,
	int* Z,               // Not presently used in kernel
	int* N,               // Not presently used in kernel
	int* Params1,
	float* Params2,
	int* NumReactingSpecies,
    int* Reactant1,
	int* Reactant2,
	int* Reactant3
) 
{	
	#define THIRD 0.333333333333333

	// Device function signatures
	__device__ int threadID (void);
	__device__ float computeTimestep (float, float, float);
	__device__ int checkAsy(float, float, float);
	__device__ float asymptoticUpdate(float, float, float, float);
	__device__ float eulerUpdate(float, float, float);
	__device__ float computekeff(float, float);
	
	// Rename integer variables passed in the array Params1
	int numberSpecies = Params1[0];
	int numberReactions = Params1[1];
	int totalFplus = Params1[2];
	int totalFminus = Params1[3];
	
	// Rename float variables passed in the array Params2
	float T9 = Params2[0];
	float tmax = Params2[1];
	float dt_init = Params2[2];

	// Logical control parameters	
	int doParallel = 1;
	int calcFluxes = 1;
	
	int intSteps = 0;  // Number of network integration steps taken
		
	// Compute the temperature-dependent factors for the rates.  Since we assume the GPU integration
	// to be done at constant temperature and density, these only need be calculated once per GPU call.
	
	float T93 = __powf(T9, THIRD);
	float t1 = 1/T9;
	float t2 = 1/T93;
	float t3 = T93;
    float t4 = T9;
    float t5 = T93*T93*T93*T93*T93;
    float t6 = __logf(T9);
	
	int i = threadID();    // Unique thread index for arbitrary number of blocks
	
	if(doParallel == 1)
	{
		// Parallel version: Compute all rates in parallel, with CUDA supplying a separate thread for 
		// each value of i.
		
		if(i < numberReactions)    // Prevent processing on threads outside bounds of network
		Rate[i] =   
			Prefac[i]* __expf(P0[i] + t1*P1[i] + t2*P2[i] + t3*P3[i] + t4*P4[i] + t5*P5[i] + t6*P6[i]); 
		
	} else {
	
		// Serial version for reference
		
		for(int i=0; i<numberReactions; i++)
		{
			Rate[i] =
			  Prefac[i]* __expf(P0[i] + t1*P1[i] + t2*P2[i] + t3*P3[i] + t4*P4[i] + t5*P5[i] + t6*P6[i]); 
		}
		
	}
	
	__syncthreads();
		
	/*
	 * Begin the time integration from t=0 to tmax. Rather than t=0 we start at some very
	 * small value of t.
	 */
	
	float t = 1.0e-16;             // The current integration time
	float dt = dt_init;            // The current integration timestep
	float prevdt = dt_init;        // The integration timestep from the previous step
		
	// Main time integration loop
	
	while(t < tmax)
	{		
		// Compute the fluxes from the previously-computed rates and the current abundances
		
		if(calcFluxes == 1)
		{
					
			if(doParallel == 1)
			{
				// Parallel version 
				
				if(i < numberReactions)

				{
					
					int nr = NumReactingSpecies[i];        // Number reacting species (1, 2, or 3)
					
					__syncthreads();

					// Switch on whether 1-body, 2-body, or 3-body reaction
					
					
					switch(nr)
					{
						
						case 1:    // 1-body; flux = rate x Y

							Flux[i] = Rate[i] * Y[*(Reactant1+i)]; 
							
						break;
							
						case 2:    // 2-body; flux = rate x Y x Y
							
							Flux[i] = Rate[i] * Y[*(Reactant1+i)] * Y[*(Reactant2+i)]; 
							
						break;
							
						case 3:    // 3-body; flux = rate x Y x Y x Y
							
							Flux[i] = Rate[i] * Y[*(Reactant1+i)] * Y[*(Reactant2+i)] 
									* Y[*(Reactant3+i)]; 
							
						break;
					}
				}
							
			} else {
		
				// Serial version for reference
			
				for(int j=0; j<numberReactions; j++)
				{
					int nr = NumReactingSpecies[j];
					
					switch(nr)
					{		
						case 1:
							Flux[j] = Rate[j] * Y[*(Reactant1+j)];			
						break;
							
						case 2:					
							Flux[j] = Rate[j] * Y[*(Reactant1+j)] * Y[*(Reactant2+j)]; 				
						break;
							
						case 3:					
							Flux[j] = Rate[j] * Y[*(Reactant1+j)] * Y[*(Reactant2+j)] 
										* Y[*(Reactant3+j)];				
						break;
					}		
				}	
			} 
		}
		
		__syncthreads();
		
		// Populate the F+ and F- arrays in parallel from the master Flux array
		
		if(i < totalFplus) 
		{
			int indy = MapFplus[i];
			Fplus[i] = FplusFac[i]*Flux[indy];
		}
		__syncthreads();
		
		 		
		if(i < totalFminus) 
		{
			Fminus[i] = FminusFac[i]*Flux[MapFminus[i]];
		}
		__syncthreads();
		
		
		// Sum the F+ and F- for each isotope.  The outer loop (in i) is parallel
		// but the inner loops (in j) are serial.  Can we do better? Seems we should
		// be able to?
			
		if(i < numberSpecies)
		{		
			// Partially serial Sum F+
			int minny = 0;
			if(i>0) minny = FplusMax[i-1]+1;
			FplusSum[i] = 0.0f;	
			for(int j=minny; j<=FplusMax[i]; j++)
			{
				FplusSum[i] += Fplus[j];
			}
				
			// Partially serial Sum F-
			minny = 0;
			if(i>0) minny = FminusMax[i-1]+1;
			FminusSum[i] = 0.0f;
			for(int j=minny; j<=FminusMax[i]; j++)
			{
				FminusSum[i] += Fminus[j];
			}
		}

		__syncthreads();
		
			
		/*
		Now use the fluxes to update the populations in parallel for this timestep
		For now we shall assume the asymptotic method. We determine whether each isotope 
		satisfies the asymptotic condition. If it does we update with the asymptotic formula. 
		If not, we update numerically using the forward Euler formula.
		*/
		
		if(i < numberSpecies)
		{		
			if(checkAsy(Fminus[i], Y[i], dt) == 1)
			{
				Y[i] = asymptoticUpdate(FplusSum[i], FminusSum[i], Y[i], dt);
			}
			else
			{
				Y[i] += eulerUpdate(FplusSum[i], FminusSum[i], dt);
			}		
		}
		__syncthreads();
		
		// Increment the integration time and set the new timestep	
		
		t += dt;
		intSteps ++;
		prevdt = dt;
		dt = computeTimestep(prevdt, t, tmax);
		
		// Temporary diagnostic halt
		if(intSteps >= 99) break;
		
	}      //+++ End of time integration while-loop +++//
	
}      // +++ End of kernel integrateNetwork +++ //


// Function to determine whether an isotope specified by speciesIndex satisfies the
// asymptotic condition. Returns 1 if it does and 0 if not.

__device__ int checkAsy(float Fminus, float Y, float dt)
{
	if(Y>0.0f && Fminus*dt/Y > 1.0f)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


// Function to return the updated Y using the asymptotic formula

__device__ float asymptoticUpdate(float Fplus, float Fminus, float Y, float dt)
{
	return (Y + Fplus*dt)/(1.0f + Fminus*dt/Y);  // Sophia He formula
}


// Function to return the Y specified by speciesIndex updated using the forward Euler method

__device__ float eulerUpdate(float FplusSum, float FminusSum, float dt)
{
	return (FplusSum-FminusSum)*dt;
}


/*
 * Construct a unique thread ID from the built-in CUDA variables threadIdx.x, blockIdx.x, and blockDim.x.
 * threadIdx.x is a unique thread ID within a block, blockIdx.x identifies the block, and blockDim.x 
 * is the number of threads in each block. (We also have available from CUDA gridDim.x, which is the 
 * number of threads in each block, but don't need it here.)  Thus the following definition ensures that the 
 * ID i is unique for every thread distributed on the device.  Note:  if all threads were assigned to a
 * single block, i=threadIdx.x would be a unique identifier.  But on a CUDA 1.1 device the maximum
 * number of threads per block is 512.  Thus, if the task requires more than 512 threads it must be
 * assigned to more than one block and i=threadIdx.x is no longer unique.
*/
	
__device__ int threadID()
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	return i;      // Return unique thread ID
}

// Function to compute the network timestep. For now it is a placeholder
// just returning the timestep as a fixed fraction of the time.

__device__ float computeTimestep(float prevdt, float t, float tmax)
{	
	float dt;
	if(t == 0.0f)
	{
		dt = 1.0e-20;
	}
	else
	{
		dt = 0.1f*t;
	}
	
	// Prevent final integration step from overstepping tmax
	//if(t+dt > tmax) dt = tmax - t;
	
	return dt;
}

// Function to compute the effective decay constant keff for asymptotic approximation
// *** NOT PRESENTLY USED ***

__device__ float computekeff(float Fminus, float Y)
{
	if(Y > 0)
	{
		return Fminus/Y;
	}
	else
	{
		return 0.0f;
	}
}



//---------- Code below is executed on the host (CPU) -----------//

#include <stdio.h>
#include<time.h>

// SIZE defines the number of reactions to be calculated.  Need min of SIZE=4395 for 365-element network, 
// 1603 for 150-isotope network, 2234 for the 194-isotope network, 3176 for the 268-isotope network,
// 48 for the alpha network, and 1566 for the nova134 network. These sizes are hardwired for now but 
// eventually we want to assign them dynamically.

#define SIZE 48                         // Max number of reactions
#define ISOTOPES 16                      // Max isotopes in network (16 for alpha network)
#define LABELSIZE 35                      // Max size of reaction string a+b>c in characters
#define PF 24                             // # entries in partition function table for each isotope
#define THIRD 0.333333333333333


// Define some CPU timing utilities. Tends to return zero for anything that takes
// less than about 10 ms.  Usage:
//     START_CPU;
//     ... code to be timed ...
//     STOP_CPU;
//     PRINT_CPU
// in the host code. These may be used to time CPU processes in the host code.

clock_t startCPU, stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("\nTimer: %f ms used by CPU\n",1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));


// Define some GPU timing utilities. These are invoked from the host program. Usage:
//     START_GPU;
//         kernelFunction <<< numBlocks, threadsPerBlock >>> (args)
//     STOP_GPU;
//     PRINT_GPU
// in the host code. This estimates the time for the kernel kernelFunction to run on the GPU.
// For a more extensive discusion, see Section 5.1.2 of the CUDA Best Practices Guide at
// http://developer.download.nvidia.com/compute/DevZone/docs/html/C/doc/CUDA_C_Best_Practices_Guide.pdf

float timeGPU;
cudaEvent_t start, stop;
#define START_GPU cudaEventCreate(&start); cudaEventCreate(&stop); cudaEventRecord(start, 0);
#define STOP_GPU cudaEventRecord(stop, 0); cudaEventSynchronize(stop);\
   cudaEventElapsedTime(&timeGPU, start, stop);\
   cudaEventDestroy(start);cudaEventDestroy(stop);
#define PRINT_GPU printf("\n\nTime to compute on GPU: %f ms \n", timeGPU);
   
// Define a utility to check for CUDA errors.  Place it immediately after a CUDA kernel
// call in the host code. The initial cudaDeviceSynchronize() command ensures that the device
// has completed all preceding requested tasks.

#define CUDA_ERROR_CHECK cudaDeviceSynchronize(); cudaError_t error = cudaGetLastError();\
   if(error != cudaSuccess){printf("***CUDA error: %s\n", cudaGetErrorString(error)); exit(-1);}\
   else{printf("\nNo CUDA errors detected\n" );}

FILE *fr;            // File pointer for data read-in

// Variables to hold data read in

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
int Reactant1[SIZE];          // Isotope index of first reactant
int Reactant2[SIZE];          // Isotope index of second reactant (if 2-body or 3-body)
int Reactant3[SIZE];		 // Isotope index of third reactant (if 3-body)

int Z[ISOTOPES];                // Array holding Z values for isotopes
int N[ISOTOPES];                // Array holding N values for isotopes
float AA[ISOTOPES];             // Array holding A values for isotopes
float Y[ISOTOPES];              // Array holding abundances Y for isotopes
float Diagnose[DIAGNOSE_SIZE];  // Diagnostics array
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

float Dens[3];

// Array with entries +1 if a reaction increases the population of the isotope (contributes to 
// F+), -1 if it decreases it (contributes to F-) and 0 if the reaction does not change the population
// of the isotope. This array is populated in the function parseF().  It is characteristic of the
// structure of the network and thus has to be calculated only once for a given network.

int reacMask[ISOTOPES][SIZE];

// Function Signatures:
void devcheck(int);
void readLibraryParams(char *);
void readNetwork(char *);
void writeNetwork(void);
void testTimerCPU(void);
void computeRatesCPU(float);
void writeRates(char *);
void writeAbundances(void);
void parseF(void);

// Filename for input rates library data. The file rateLibrary.data output by the Java code through 
// the stream toRateData has the expected format for this file. Standard test cases: 
// rateLibrary_alpha.data, rateLibrary_150.data, rateLibrary_365.data, rateLibrary_nova134.data.

char rateLibraryFile[] = "rateLibrary_alpha.data";

// Filename for network + partition function input.  The file output/CUDAnet.inp
// output by the Java code through the stream toCUDAnet has the expected format for 
// this file. Standard test cases: CUDAnet_alphasolar.inp, CUDAnet_150solar.inp,
// CUDAnet_365solar.inp, CUDAnet_nova134.inp.

char networkFile[] = "CUDAnet_alpha.inp";

// Control some of the printout (1 to include, 0 to suppress)
int displayInput = 0;
int ratesCPU = 0;

// Total number of F+ and F- terms in the network
int totalFplus = 0;
int totalFminus = 0;

// Arrays to hold non-zero fluxes in the network. These will be allocated dynamically 
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
// Fminus arrays. These will be allocated dynamically below malloc

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

// Control diagnostic printout (1 to print, 0 to suppress)
int showParsing = 0;
int showFparsing = 0;

// Array of integer input parameters.  Needed because we can't pass
// too many arguments to the kernel (256 byte limit for 1.1 devices).

int Params1[4];

// Array of float input parameters.  Needed because we can't pass
// too many arguments to the kernel (256 byte limit for 1.1 devices).

float Params2[3];


// Main CPU routine

int main()
{
	// Ensure that a valid device (GPU) exists 
	printf("\nChecking for valid device:\n");
	devcheck(0);
	
	// Check available memory on the GPU	
	size_t msizeFree;
	size_t msizeTotal;	
	cudaMemGetInfo(&msizeFree, &msizeTotal);	 
	printf("\nGPU total memory: %d\nGPU free memory: %d", (int)msizeTotal, (int)msizeFree);
	
	// Following memory queries not supported on GF 8600 GT	
	//  size_t msize;
	// 	cudaDeviceGetLimit(&msize, cudaLimitMallocHeapSize);
	//  printf("\nGPU heap size: %d", (int)msize);
	// 	cudaDeviceGetLimit(&msize, cudaLimitStackSize);
	// 	printf("\nGPU stack size for each thread: %d\n", (int)msize);
	
	// Set the temperature in units of 10^9 K and density in units of g/cm^3. The 
	// temperature and density will be passed from the hydro code in an operator-split 
	// coupling of this network to hydro. These will be used to calculate the reaction
	// rates in the network on the GPU. Since we are assuming operator splitting, the
	// temperature and density are assumed constant for the entire network integration
	// on the gPU.
	
	float T9 = 6.0f;
	float rho = 1.0e8;
	
	// Set the range of time integration and the initial timestep.  In an operator-split
	// coupling tmax will come from the hydro and dt_init will likely be the last timestep
	// of the previous network integration (for the preceding hydro timestep).
	
	float tmax = 1e-11;
	float dt_init = 1e-17;               
	
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
	
	// Optionally compute the temperature-dependent rates first on the CPU
	if(ratesCPU == 1){
		
		// Test the CPU timer by executing a long, pointless loop
		testTimerCPU();
		
		START_CPU     // Start a timer for the actual calculation
		
		// First compute the rates serially for reference using the CPU
			
		computeRatesCPU(T9);
		
		STOP_CPU;     // Stop the timer
		PRINT_CPU;    // Print timing information
		
		// Display the rates. (Note the two different techniques used in calling
		// writeRates here and for "GPU" below). Here we use a pointer; in the
		// example below we pass an array.)
		
		char label[] = "on CPU";
		char *labelPtr = label;
		writeRates(labelPtr);
	
	}
	
	// Read in network file and associated partition functions.  This is required only
	// once at the very beginning of the calculation.
	
	char *networkFilePtr = networkFile;
	readNetwork(networkFilePtr);
	writeNetwork();

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
	
	// Create 1D arrays to hold non-zero F+ and F- for all reactions for all isotopes,
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
	// Fplus and Fminus 1D arrays. Since FplusMax and FplusMin are related, and likewise
	// FminusMax and FminusMin are related, we will only need to pass FplusMax and
	// FminusMax to the kernel.
	
	FplusMax = (int*) malloc(sizeof(int) * numberSpecies);
	FplusMin = (int*) malloc(sizeof(int) * numberSpecies);
	FminusMax = (int*) malloc(sizeof(int) * numberSpecies);
	FminusMin = (int*) malloc(sizeof(int) * numberSpecies);
	
	// Fill the integer parameter array to pass to the kernel. Doing this to bypass
	// the formal limit of 256 for arguments passed to the kernel.
	
	Params1[0] = numberSpecies;
	Params1[1] = numberReactions;
	Params1[2] = totalFplus;
	Params1[3] = totalFminus;
	
	// Fill the float parameter array to pass to the kernel. Doing this to bypass
	// the formal limit of 256 for arguments passed to the kernel.
	
	Params2[0] = T9;
	Params2[1] = tmax;
	Params2[2] = dt_init; 		

	// Create 1D arrays that will be used to map finite F+ and F- to the Flux array.
	
	FplusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
	FminusIsotopeCut = (int*) malloc(sizeof(int) * numberSpecies);
	
	FplusIsotopeIndex = (int*) malloc(sizeof(int) * totalFplus);
	FminusIsotopeIndex = (int*) malloc(sizeof(int) * totalFminus);
	
	// Create 1D arrays that will hold the index of the isotope for the F+ or F- term
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
	
	// Diagnostic output
	if(showFparsing == 1)
	{
		printf("\n\n\nMAX F+ and F- INDEX FOR EACH ISOTOPE:\n");	
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
	// number of occurences of the species in the reaction.  Note that this can only
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
	
	// Diagnostic output
	
	if(showFparsing == 1)
	{
		printf("\n\n\n---------- %d NON-VANISHING F+ SOURCE TERMS ----------\n", totalFplus);
		printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
					isoLabel[FplusIsotopeIndex[0]], FplusIsotopeIndex[0],
					numFluxPlus[FplusIsotopeIndex[0]]);
		for(int i=0; i<totalFplus; i++)
		{
			printf("\n   Isotope index = %d F+ index = %d Reac index = %d  %s", 
					FplusIsotopeIndex[i], i,
				MapFplus[i], reacLabel[MapFplus[i]]); 
			if(i == (FplusIsotopeCut[FplusIsotopeIndex[i]] - 1)  && i != totalFplus-1)
			{
				printf("\n");
				printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):", 
						isoLabel[FplusIsotopeIndex[i+1]], FplusIsotopeIndex[i+1],
						numFluxPlus[FplusIsotopeIndex[i+1]]);
			}
		}	
		
		printf("\n\n\n---------- %d NON-VANISHING F- SOURCE TERMS ----------\n", totalFminus);
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
		
		printf("\n");
	}
	
	// Set up the device pointers corresponding to the arrays.  Required only once at the
	// very beginning of the calculation.
	
	float *devPtrP0;
	float *devPtrP1; 
	float *devPtrP2; 
	float *devPtrP3;
	float *devPtrP4;
	float *devPtrP5;
	float *devPtrP6;
	float *devPtrPrefac;
	float *devPtrQ;
	float *devPtrRate;
    float *devPtrFlux;
	float *devPtrFplus;
	float *devPtrFminus;
	float *devPtrFplusFac;
	float *devPtrFminusFac;
	float *devPtrFplusSum;
	float *devPtrFminusSum;
	int *devPtrFplusMax;
	int *devPtrFminusMax;
	int *devPtrMapFplus;
	int *devPtrMapFminus;
	float *devPtrY;
	float *devPtrDiagnose;
	int *devPtrZ;                       // Z is not presently used in the kernel
	int *devPtrN;                       // N is not presently used in the kernel
	int *devPtrParams1;
	float *devPtrParams2;
	int *devPtrNumReactingSpecies;
    int *devPtrReactant1;
	int *devPtrReactant2;
	int *devPtrReactant3;
	
	// Allocate float and int memory on the device (the GPU) for all variables.  Required only
	// once at the very beginning of the calculation.
	
	int memsize1 = SIZE*sizeof(float);         // Memory size for floats labeled by reaction index
	int memsize2 = SIZE*sizeof(int);           // Memory size for ints labeled by reaction index
	int memsize3 = ISOTOPES*sizeof(float);     // Memory size for floats labeled by isotope index
	int memsize4 = ISOTOPES*sizeof(int);       // Memory size for ints labeled by isotope index
	int memsize5 = totalFplus*sizeof(float);   // Memory size for contributing F+ values
	int memsize6 = totalFminus*sizeof(float);  // Memory size for contributing F- values
	int memsize7 = totalFplus*sizeof(int);     // Memory size for contributing F+ value indices
	int memsize8 = totalFminus*sizeof(int);    // Memory size for contributing F- value indices
	
	cudaMalloc((void**)&devPtrP0, memsize1); 
	cudaMalloc((void**)&devPtrP1, memsize1); 
	cudaMalloc((void**)&devPtrP2, memsize1);
	cudaMalloc((void**)&devPtrP3, memsize1);
	cudaMalloc((void**)&devPtrP4, memsize1);
	cudaMalloc((void**)&devPtrP5, memsize1);
	cudaMalloc((void**)&devPtrP6, memsize1);
	cudaMalloc((void**)&devPtrPrefac, memsize1);
	cudaMalloc((void**)&devPtrQ, memsize1);
	cudaMalloc((void**)&devPtrRate, memsize1);
    cudaMalloc((void**)&devPtrFlux, memsize1);
	cudaMalloc((void**)&devPtrFplus, memsize5);
	cudaMalloc((void**)&devPtrFminus, memsize6);
	cudaMalloc((void**)&devPtrFplusFac, memsize5);
	cudaMalloc((void**)&devPtrFminusFac, memsize6);
	cudaMalloc((void**)&devPtrFplusSum, memsize3);
	cudaMalloc((void**)&devPtrFminusSum, memsize3);
	cudaMalloc((void**)&devPtrFplusMax, memsize4);
	cudaMalloc((void**)&devPtrFminusMax, memsize4);	
	cudaMalloc((void**)&devPtrMapFplus, memsize7);
	cudaMalloc((void**)&devPtrMapFminus, memsize8);
	cudaMalloc((void**)&devPtrY, memsize3);
	cudaMalloc((void**)&devPtrDiagnose, DIAGNOSE_SIZE*sizeof(float));
	cudaMalloc((void**)&devPtrZ, memsize4);
	cudaMalloc((void**)&devPtrN, memsize4);
	cudaMalloc((void**)&devPtrParams1, 4*sizeof(int));
	cudaMalloc((void**)&devPtrParams2, 3*sizeof(float));
	cudaMalloc((void**)&devPtrNumReactingSpecies, memsize2);
    cudaMalloc((void**)&devPtrReactant1, memsize2);
	cudaMalloc((void**)&devPtrReactant2, memsize2);
	cudaMalloc((void**)&devPtrReactant3, memsize2);
	
	/*
	 * Copy array memory to the device using cudaMemcpy(void* destination, 
	 * const void* source, size_t memSize, enum cudaMemcpyKind kind), which 
	 * copies memSize bytes from the memory area pointed to by source to the 
	 * memory area pointed to by destination, with kind=cudaMemcpyHostToDevice 
	 * or kind=cudaMemcpyDeviceToHost, or kind=cudaMemcpyDeviceToDevice 
	 * specifying the nature (direction) of the copy. Note: in C an array name
	 * is effectively a pointer.  This is required only once at the very beginning
	 * of the calculation. These quantities will then reside on the GPU for the
	 * duration of the calculation (over all hydro timesteps).
	*/
	
	cudaMemcpy(devPtrP0, P0, memsize1, cudaMemcpyHostToDevice); 
	cudaMemcpy(devPtrP1, P1, memsize1, cudaMemcpyHostToDevice); 
	cudaMemcpy(devPtrP2, P2, memsize1, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrP3, P3, memsize1, cudaMemcpyHostToDevice); 
	cudaMemcpy(devPtrP4, P4, memsize1, cudaMemcpyHostToDevice); 
	cudaMemcpy(devPtrP5, P5, memsize1, cudaMemcpyHostToDevice); 
	cudaMemcpy(devPtrP6, P6, memsize1, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrPrefac, Prefac, memsize1, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrQ, Q, memsize1, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFplus, Fplus, memsize5, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFminus, Fminus, memsize6, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFplusFac, FplusFac, memsize5, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFminusFac, FminusFac, memsize6, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFplusSum, FplusSum, memsize3, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFminusSum, FminusSum, memsize3, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFplusMax, FplusMax, memsize4, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrFminusMax, FminusMax, memsize4, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrMapFplus, MapFplus, memsize7, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrMapFminus, MapFminus, memsize8, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrY, Y, memsize3, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrDiagnose, Diagnose, DIAGNOSE_SIZE*sizeof(float), cudaMemcpyHostToDevice);  // Needed?
	cudaMemcpy(devPtrZ, Z, memsize4, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrN, N, memsize4, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrParams1, Params1, 4*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrParams2, Params2, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrNumReactingSpecies, NumReactingSpecies, memsize2, cudaMemcpyHostToDevice);
    cudaMemcpy(devPtrReactant1, Reactant1, memsize2, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrReactant2, Reactant2, memsize2, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrReactant3, Reactant3, memsize2, cudaMemcpyHostToDevice);

	/*
	 * Define the number of blocks on the grid and the number of
	 * threads per block. The number of reactions is given by SIZE, so
	 * the product of numBlocks.x and threadsPerBlock.x must be greater
	 * than or equal to SIZE. If we choose threadsPerBlock.x = 256,
	 * we require at least numBlocks.x = 1 for the alpha network, 7 for the 150-
	 * isotope network, 9 for the 194-isotope network, 13 for the 268-
	 * isotope network, 18 for the 365-isotope network, and 7 for the
	 * nova134 network. Being done by hand now, but eventually we should 
	 * automate this. Required only once at the very beginning of the calculation.
	*/
	
	
	dim3 numBlocks(4, 1, 1);
	dim3 threadsPerBlock(BLOCKSIZE, 1, 1);
	
	/* 
	 * Launch kernel integrateNetwork on GPU. This will execute the global function 
	 * integrateNetwork on the device using numBlocks blocks, each with threadsPerBlock 
	 * threads.  Pass pointers for the arrays and values for the scalars.  This kernel 
	 * will be launched once for every hydro timestep, with the entire network calculation 
	 * for the hydro timestep being done on the GPU within this kernel. At the beginning 
	 * we must pass to the kernel running the network on the GPU the current temperature, 
	 * the proposed initial network timestep (likely we will choose the timestep from the
	 * end of the network integration for the preceding hydro timestep), and the current 
	 * abundances Y[i] for the species in the network (which may have been altered in the 
	 * hydro timestep by processes like advection. 
	*/
	
	START_GPU;     // Start timer for device code
	
	integrateNetwork <<< numBlocks, threadsPerBlock >>> 
	(
		devPtrP0, 
		devPtrP1, 
		devPtrP2, 
		devPtrP3, 
		devPtrP4, 
		devPtrP5, 
		devPtrP6, 
		devPtrPrefac,  
		devPtrQ, 
		devPtrRate, 
		devPtrFlux, 
		devPtrFplus,
        devPtrFminus,
		devPtrFplusFac,
		devPtrFminusFac,
		devPtrFplusSum,
		devPtrFminusSum,
 		devPtrFplusMax,
 		devPtrFminusMax,
        devPtrMapFplus,
		devPtrMapFminus,
		devPtrY, 
		devPtrDiagnose,
		devPtrZ, 
		devPtrN,
		devPtrParams1,
		devPtrParams2,
        devPtrNumReactingSpecies,
		devPtrReactant1, 
		devPtrReactant2, 
		devPtrReactant3
	);

	STOP_GPU;           // Stop timer for device code
	PRINT_GPU;          // Print timing for device code
	
	CUDA_ERROR_CHECK    // Check for CUDA errors
	
	// The integration is finished on the GPU.  Copy the results from the device back 
	// to the host.
	
	cudaMemcpy(Rate, devPtrRate, memsize1, cudaMemcpyDeviceToHost);
    cudaMemcpy(Flux, devPtrFlux, memsize1, cudaMemcpyDeviceToHost);
	cudaMemcpy(Fplus, devPtrFplus, memsize5, cudaMemcpyDeviceToHost);
	cudaMemcpy(Fminus, devPtrFminus, memsize6, cudaMemcpyDeviceToHost);
	cudaMemcpy(FplusFac, devPtrFplusFac, memsize5, cudaMemcpyDeviceToHost);
	cudaMemcpy(FminusFac, devPtrFminusFac, memsize6, cudaMemcpyDeviceToHost);
	cudaMemcpy(FplusSum, devPtrFplusSum, memsize3, cudaMemcpyDeviceToHost);
	cudaMemcpy(FminusSum, devPtrFminusSum, memsize3, cudaMemcpyDeviceToHost);
	cudaMemcpy(FplusMax, devPtrFplusMax, memsize4, cudaMemcpyDeviceToHost);
	cudaMemcpy(FminusMax, devPtrFminusMax, memsize4, cudaMemcpyDeviceToHost);
	cudaMemcpy(MapFplus, devPtrMapFplus, memsize7, cudaMemcpyDeviceToHost);
	cudaMemcpy(MapFminus, devPtrMapFminus, memsize8, cudaMemcpyDeviceToHost);
	cudaMemcpy(Y, devPtrY, memsize3, cudaMemcpyDeviceToHost);
	cudaMemcpy(Diagnose, devPtrDiagnose, DIAGNOSE_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	
	// Display the results returned from the GPU
    char label2[] = "on GPU";
    writeRates(label2);
	
	// Print the final values of F+ and F- transferred from GPU
	printf("\n\nFINAL F+ VALUES:\n");
	for(int i=0; i<totalFplus; i++)
	{
		printf("\nF+[%d] = %7.4e  Increases Y[%s] through %s  MapIndex=%d  FplusFac=%3.1f", 
			i, Fplus[i], isoLabel[FplusIsotopeIndex[i]], 
			reacLabel[MapFplus[i]], MapFplus[i], FplusFac[i]);
	}
	printf("\n\n\nFINAL F- VALUES:\n");
	for(int i=0; i<totalFminus; i++)
	{
		printf("\nF-[%d] = %7.4e  Decreases Y[%s] through %s  MapIndex=%d  FminusFac=%3.1f", 
			i, Fminus[i], isoLabel[FminusIsotopeIndex[i]], 
			reacLabel[MapFminus[i]], MapFminus[i], FminusFac[i]);
	}
		
	printf("\n\n\nF+ and F- MIN AND MAX FOR EACH ISOTOPE:\n");
	for(int i=0; i<numberSpecies; i++)
	{
		printf("\n%3d %5s F+min=%3d F+max=%3d F-min=%3d F-max=%d", 
			   i, isoLabel[i], FplusMin[i], FplusMax[i], FminusMin[i], FminusMax[i]);
	}
	
	printf("\n\n\nSUM OF FLUXES FOR EACH ISOTOPE:\n");
	
	// Print the total F+ and F- for each isotope transferred from the GPU
	float totalFplus = 0.0f;
	float totalFminus = 0.0f;
	for(int i=0; i<numberSpecies; i++)
	{
		printf("\n%3d %5s  sumF+=%10.4e  sumF-=%10.4e Fnet=%10.4e Y=%10.4e", 
			   i, isoLabel[i], FplusSum[i], FminusSum[i], FplusSum[i]-FminusSum[i], Y[i]);
		totalFplus += FplusSum[i];
		totalFminus += FminusSum[i];
	}
	
	printf("\n\ntotalF+ = %7.4e  totalF- = %7.4e", totalFplus, totalFminus);

	// Diagnostics returned from the GPU
	
	printf("\n\n\nDIAGNOSTICS:\n\n");
	for(int i=0; i<DIAGNOSE_SIZE; i++)
	{
		printf("Diagnose[%d]=%7.4e\n", i, Diagnose[i]);
	}
	
	printf("\n\nFINAL ABUNDANCES:\n");
	writeAbundances();
	
	PRINT_GPU;

	printf("\n");
	
    // Free memory allocated on the device  
	
	cudaFree(devPtrP0); 
	cudaFree(devPtrP1); 
	cudaFree(devPtrP2);
	cudaFree(devPtrP3);
	cudaFree(devPtrP4);
	cudaFree(devPtrP5);
	cudaFree(devPtrP6);
	cudaFree(devPtrPrefac);
	cudaFree(devPtrQ);
	cudaFree(devPtrRate);
    cudaFree(devPtrFlux);
	cudaFree(devPtrFplus);
	cudaFree(devPtrFminus);
	cudaFree(devPtrFplusFac);
	cudaFree(devPtrFminusFac);
	cudaFree(devPtrFplusSum);
	cudaFree(devPtrFminusSum);
	cudaFree(devPtrFplusMax);
	cudaFree(devPtrFminusMax);
	cudaFree(devPtrMapFplus);
	cudaFree(devPtrMapFminus);
	cudaFree(devPtrY);
	cudaFree(devPtrDiagnose);
	cudaFree(devPtrZ);
	cudaFree(devPtrN);
	cudaFree(devPtrParams1);
	cudaFree(devPtrParams2);
	cudaFree(devPtrNumReactingSpecies);
    cudaFree(devPtrReactant1);
	cudaFree(devPtrReactant2);
	cudaFree(devPtrReactant3);
	
	// Free memory allocated on the CPU
	
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
	
}  // End main


// Function to compute the rates on the CPU rather than GPU.

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
}

// Function to print out all the rates. The label can be used to distinguish cases
// if called more than once.

void  writeRates(char *label)
{
	printf("\nCOMPUTED RATES (%s):\n\n", label);
	for (int i=0; i<numberReactions; i++) 
	{
		printf("%d %s Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Flux=%7.4e Q=%6.3f Prefac=%6.3e Reactants=%d\n",
			i,reacLabel[i], Rate[i], Y[Reactant1[i]], Y[Reactant2[i]], Y[Reactant3[i]], Flux[i], 
			Q[i], Prefac[i], NumReactingSpecies[i]);
	}
}

// Function to write out the abundances in the network

void writeAbundances()
{
	printf("\nIndex  Isotope   Abundance Y   Mass Frac X");
	float sumX = 0.0f;
	for(int i=0; i<ISOTOPES; i++)
	{
		float X = Y[i]*AA[i];
		sumX += X;
		printf("\n %4d     %4s    %8.4e    %8.4e", i, isoLabel[i], Y[i], X);
	}
	printf("\n\nsum X = %6.4f", sumX);
}

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
	printf("\n T9 = %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
		Tpf[0],Tpf[1],Tpf[2],Tpf[3],Tpf[4],Tpf[5],Tpf[6],Tpf[7],
		Tpf[8],Tpf[9],Tpf[10],Tpf[11],Tpf[12],Tpf[13],Tpf[14],Tpf[15],
		Tpf[16],Tpf[17],Tpf[18],Tpf[19],Tpf[20],Tpf[21],Tpf[22],Tpf[23]
	);
	for(int j=0; j<numberSpecies; j++){
		printf("\n%-5s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
	isoLabel[j],pf[j][0],pf[j][1],pf[j][2],pf[j][3],pf[j][4],pf[j][5],pf[j][6],pf[j][7],	
	pf[j][8],pf[j][9],pf[j][10],pf[j][11],pf[j][12],pf[j][13],pf[j][14],pf[j][15],
	pf[j][16],pf[j][17],pf[j][18],pf[j][19],pf[j][20],pf[j][21],pf[j][22],pf[j][23] );
	}
	
 	printf("\n");
}


/* Function to read rate parameter data file line by line, with filename as argument.
 This file is expected to have one reaction per line with the line structure
	 p0 p1 p2 p3 p4 p5 p6 reactionLabel
 where the pn are the values of the 7 Reaclib parameters for a reaction,
 reactionLabel is a label for the reaction that must contain no whitespace, and
 all fields on a line are separated by a blank space.
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
	Read in the file line by line and parse into variables.  The expected
	structure of each line is
	     float float float float float float float string
	each separated by a space, with no whitespace in the string.
	(See http://stackoverflow.com/questions/2854488/reading-a-string-with-spaces-with-sscanf
	for how to read string with spaces.)
	*/
	
	int n = -1;
	int subindex = -1;
	
	if(displayInput == 1) printf("\nData read in:\n\n");
	
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
	
	for(int i=0; i<numberReactions; i++)
	{
		Reactant1[i] = ReactantIndex[i][0];
		Reactant2[i] = ReactantIndex[i][1];
		Reactant3[i] = ReactantIndex[i][2];
	}
	
	fclose(fr);           // Close the file
}


/* Function to read the network data file line by line, with the filename as argument.
 This file is expected to have 4 lines per isotope with the line structure
	 isotopeSymbol A  Z  N  Y  MassExcess
	 pf00 pf01 pf02 pf03 pf04 pf05 pf06 pf07
	 pf10 pf11 pf12 pf13 pf14 pf15 pf16 pf17
	 pf20 pf21 pf22 pf23 pf24 pf25 pf26 pf27
where isotopeSymbol is an isotope label, A=Z+N is the atomic mass number, Z is the proton number, 
N is the neutron number, Y is the current abundance, MassExcess is the mass
excess in MeV, and the pf are 24 values of the partition function for that isotope at
different values of the temperature that will form a table for interpolation in temperature.
The assumed 24 values of the temperature for the partition function table are
{ 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 
4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 } in units of 10^9 K.
All fields on a line are separated by a blank space and there is no whitespace in the isotopeSymbol.
The type signature of these four lines corresponding to a single isotope is
	string int int int float float
	float float float float float float float float
	float float float float float float float float
	float float float float float float float float
Here is an example for two isotopes:

ca40 40 20 20 0.0 -34.846
1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.01 1.04 1.09 1.2 1.38
ti44 44 22 22 0.0 -37.548
1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.01 1.03 1.08 1.14
1.23 1.35 1.49 1.85 2.35 3.01 3.86 4.94

A file with this format is written from the Java code to the file output/CUDAnetwork.inp using the
Java stream toCUDAnet.
	
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
	
	if(displayInput==1) printf("\nData read in:\n");
	
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
}


/*
  Function to find the contributions to F+ and F- of each reaction for each isotope.  
  This is executed only once at the beginning of the entire calculation to determine 
  the structure of the network.
*/

void parseF()
{
	if(showParsing == 1)
		printf("\nUse parseF() to find F+ and F- flux components for each species:");
	
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
	
	printf("\n\nPART OF FLUX-ISOTOPE COMPONENT ARRAY (-n --> F-; +n --> F+ for given isotope):");
	
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
	}
	
	printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out of %d x %d = %d possibilities.", 
		totalFplus, totalFminus, SIZE, ISOTOPES, SIZE*ISOTOPES);
}


// Function to check that a valid device exists. Copied from
// http://www.ncsa.illinois.edu/UserInfo/Training/Workshops/CUDA/presentations/tutorial-CUDA.html

void devcheck(int gpudevice) 
{ 
	int device_count=0; 
	int device;  // used with  cudaGetDevice() to verify cudaSetDevice() 

	// get the number of non-emulation devices  detected 
	cudaGetDeviceCount( &device_count); 
	if (gpudevice > device_count) 
	{ 
		printf("gpudevice >=  device_count ... exiting\n"); 
		exit(1); 
	} 
	cudaError_t cudareturn; 
	cudaDeviceProp deviceProp; 
    
	// cudaGetDeviceProperties() is also  demonstrated in the deviceQuery/ example
	// of the sdk projects directory 
	
	cudaGetDeviceProperties(&deviceProp,  gpudevice); 
	printf("\n[deviceProp.major.deviceProp.minor] = [%d.%d]\n", 
	deviceProp.major, deviceProp.minor); 

	if (deviceProp.major > 999) 
	{ 
		printf("warning, CUDA Device  Emulation (CPU) detected, exiting\n"); 
		exit(1); 
	} 
   
	// choose a cuda device for kernel  execution 
	cudareturn=cudaSetDevice(gpudevice); 
	if (cudareturn == cudaErrorInvalidDevice) 
	{ 
		perror("cudaSetDevice returned  cudaErrorInvalidDevice"); 
	} 
	else 
	{ 
		// double check that device was properly selected 
		cudaGetDevice(&device); 
		printf("cudaGetDevice()=%d\n",device); 
	} 
}
