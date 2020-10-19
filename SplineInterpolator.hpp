#include <cstdio>
/*
Class to implement 1D and 2D cubic spline interpolation. Adapted from algorithms 
in Numerical Recipes. For 1D interpolations, use the method spline to set up an 
interpolation table and then use the method splint to interpolate in the independent 
variable. For 2D interpolations, use the method spline2 to set up the 2D interpolation 
table and then use the method splint2 to interpolate using that table. The class also 
makes available the utility method bisection, which finds the indices in a 1-D table 
that bracket a particular entry in the table (assuming the entries to increase 
monotonically).  See the class SplineTester for examples of using SplineInterpolator.
*/


// --------------------------------------------------------------------------
//  Class SplineInterpolator implements cubic spline interpolation in one 
//  dimension and bicubic spline interpolation in two dimensions.
// --------------------------------------------------------------------------

class SplineInterpolator {

	public:

    double *x;
	double *x1a;
	double *x2a;
	double **ya;
	double *y;
	double *y2;
	double **y2a;
	double *u;
	int n,m;
	double maxx1,minx1,maxx2,minx2;
	
	SplineInterpolator() { }


	/*------------------------------------------------------------------------------
     Adaptation of 1D cubic spline algorithm from Numerical Recipes in Fortran.
     This method processes the array to compute and store the second derivatives
     that will be used to do later spline interpolation.  It assumes
     the existence of an array xarray of independent variables and an array 
	 yarray(x) of dependent variables, with the xarray array monotonically 
	 increasing.  The second derivatives are stored in the array y2.  Elements of
	 the y2 array can be accessed for diagnostic purposes using the method 
	 getSplined2(index).
    -------------------------------------------------------------------------------*/
	
	void spline(double xarray[], double yarray[]) {
		int n = arrayLength1d(xarray);
		int m = arrayLength1d(yarray);

		if(n != m){
			printf("\nWarning: lengths of x and y(x) arrays not equal: xarray.length= %d, yarray.length= %d", n, m);
		}

		//x[n] = {};
		//y[n] = {};
		//y2[n] = {};
		//u[n] = {};
		
		// Copy passed arrays to internal arrays

		for(int i=0; i<n; i++){
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
			u[i] = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])
					/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		}

		y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);

		// Backsubstitution loop of tridiagonal algorithm

		for (int i = n-2; i >= 0; i--){
			y2[i] = y2[i]*y2[i+1] + u[i];
		}
	}


	/*-------------------------------------------------------------------------
     Method splint calculates the 1D cubic spline polynomial for an arbitrary
	 argument xvalue once the second derivative table has been constructed
	 using the method spline.  This method returns the interpolated value 
	 y = f(xvalue) and may be called any number of times once the second
	 derivative table has been created.
 	---------------------------------------------------------------------------*/

	double splint(double xvalue) {

		int n = arrayLength1d(x);

		// Return -1 with error message if argument out of table bounds

		if (xvalue < x[0] || xvalue > x[n-1]) {
			printf("Argument (%f) Out of Table Bounds %f - %f", x[0], x[n-1], xvalue);
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
	 For an array xarray[] and argument xvalue, public method bisection finds 
	 the indices ilow and ihigh=ilow+1 that bracket the position of xvalue 
	 in the array.  (Method returns the value of ilow, from which 
	 ihigh = ilow+1).  The array is assumed to be monotonically increasing
	 in value. Sample method call:

			double [] x = {10,12,14,16,18};
			double xvalue = 12.2;
			int ilow = instance.bisection(x,xvalue);
			int ihigh= ilow + 1;

  	 In this case, the method returns ilow = 1 and therefore ihigh = 2.  The method
	 checks first that xvalue is within the bounds of the table and returns -1 if
	 this condition is not satisfied.
	-------------------------------------------------------------------------------*/


	int bisection(double xarray[], int arrLen, double xvalue){

		int n = arrLen;

		// Check that xvalue is within bounds of the table.  If not, quit
		// with error message and return -1

		double minx = xarray[0];
		double maxx = xarray[n-1];
		if(xvalue > maxx || xvalue < minx){
			printf("Abort bisection: argument (%f) Out of Bounds", xvalue);
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


	
	/*-----------------------------------------------------------------------------
	 Method to set up spline 2nd deriv table for a 2D spline interpolation.
	 x1array and x2array are 1D arrays holding the independent (x1,x2) variables. 
	 yarray is a 2D array holding y=f(x1,x2), however 2D arrays of undefined size
	 cannot be directly passed and so we will use a double pointer that holds the
	 address of the beginning of the 2D array. The calling function will also need 
	 to utilize double pointers to call this function. Note that the independent variable 
	 arrays are assumed to both be in ascending order.  The second derivatives are 
	 stored in the array y2a.  Elements of the y2a array can be accessed for 
	 diagnostic purposes using the method getSpline2d2(row,column).
	------------------------------------------------------------------------------*/
	
	void spline2( double x1array[], int x1Len, double x2array[], int x2Len, double **yarray) {

		n = arrayLength1d(x1array);
		m = arrayLength1d(x2array);

		//  Store min and max for later checks

		maxx1 = x1array[n-1];  
		minx1 = x1array[0];
		maxx2 = x2array[m-1];
		minx2 = x2array[0];

		// Set up temporary 1D array required for 2D interpolations.
		
		double ytmp[m]= {};
		
		double x1a[n] = {};     // 2D independent variable x1
		double x2a[m] = {};     // 2D independent variable x2

		// Set up the 2d array ya[i][j] that will hold the dependent 
		// variable y(x1,x2)

		ya = new double*[n];
		for(int i=0; i<n; i++){
			ya[i] = new double[m];
		}

		// Set up the 2d array y2a[i][j] that will hold the second derivatives

		y2a = new double*[n];
		for(int i=0; i<n; i++){
			y2a[i] = new double[m];
		}

		// Fill the local arrays x1a, x2a, and ya with the user-supplied dependent 
		// and independent variable data from arrays on the timeline containing this clip

		for (int i=0; i<n; i++) { x1a[i] = x1array[i]; }  // x1
		for (int j=0; j<m; j++) { x2a[j] = x2array[j];}   // x2
		for (int i=0; i<n; i++) {
			for(int j=0; j<m; j++) { ya[i][j] = yarray[i][j];}  // y(x1,x2)
		}

		// Compute the 2nd derivative array required for 2D interpolation

		for (int j=0; j<m; j++) {
			for (int k=0; k<n; k++) {
				ytmp[k] = ya[j][k];
			}
			// Call 1D spline; returns 2nd derivs in y2
			spline(x2a, ytmp);    
			//  Fill 2nd derivative array
			for (int k=0; k<n; k++) {
				y2a[j][k] = y2[k];        
			}
		}	
	}
	

	/*------------------------------------------------------------------------
     Method splint2 calculates the 2D cubic spline polynomial for a arbitrary
	 arguments (x1,x2) once the second derivative table has been constructed
	 using the method spline2.  This function returns the interpolated value 
	 y = f(x1,x2) and may be called any number of times once the second
	 derivative table has been created.
 	-------------------------------------------------------------------------*/

	double splint2 (double x1, double x2) {

		//  Check that the arguments lie within the range of the table.  If not,
		//  terminate with error message and return -1.

		if (x1<minx1 || x1 > maxx1){
			printf("*** Error: (x1=%f out of table range; return -1) ***",x1);
			return -1; 
		}
		if (x2<minx2 || x2 > maxx2){
			printf("*** Error: (x2=%f out of table range; return -1) ***",x2);
			return -1; 
		}
		
		// Set up temporary 1D array that will be required in 2D interpolations. 
		
		double yytmp[m] = {};

		for (int j=0; j<m; j++) {
			for (int k=0; k<n; k++) {
				y[k] = ya[j][k];
				y2[k] = y2a[j][k];
			}
			yytmp[j] = splint(x2);
		}
		spline(x1a, yytmp);
		return splint(x1);
	}


	// Public method to return element of spline 2nd derivative array for 1D
	// interpolation.  Mostly useful for diagnostics.

	double getSplined2(int index){
		return y2[index];
	}


    // Public method to return element of spline 2nd derivative array for 2D
	// interpolation.  Mostly useful for diagnostics.

	double getSpline2d2(int row, int column){
		return y2a[row][column];
	}

	int arrayLength1d(double arr[]){
		return sizeof(arr)/sizeof(arr[0]);
		//return *(&arr+1) - arr;
	}
};
