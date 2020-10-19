#include "SplineInterpolator.hpp"

// ---------------------------------------------------------------------------
//  Class SplineTester to test spline algorithms in class SplineInterpolator
// ---------------------------------------------------------------------------

class SplineTester {

	public:

    void test() {

        SplineInterpolator si = SplineInterpolator();

		// Test bisection method

		double x[] = {10,12,14,16,18};
		double xvalue = 12.2;
		int ilow = si.bisection(x, 5, xvalue);
		int ihigh= ilow + 1;

		if(ilow == -1){
			printf("\nBisection test:  Out of range; -1 returned\n");
		} else {
			printf("Bisection test:  ilow=%d ihigh=%d\n\n", ilow, ihigh);
		}

		// Test 1D spline interpolation

		double xx[] = {10,12,14,16,18};
		double yy[] = {1.3,1.8,2.0,1.6,-1.2};

		si.spline(xx,yy);

		double tryx = 12.2;
		double value = si.splint(tryx);
		printf("1D interpolation test:  y(x) = y(%f)=%f\n", tryx, value);


		// Test 2D spline interpolation

		//  Define the two 1D independent variable arrays.  In this example,
		//  just set them to one-half the row and column indices

		int rows = 11;
		int columns = 11;
		double x1[rows];
		double x2[columns];

		for (int i=0; i<rows; i++) {
			x1[i] = (double)i/2;
		}
		for (int j=0; j<columns; j++) {
			x2[j] = (double)j/2;
		}

		// Create a 2d array to hold data

		//double array2d[rows][columns];
		double **array2d = new double*[rows];

		//  Print out the independent variable arrays

		for(int i=0; i<rows; i++){
			printf("2D x1=%f", x1[i]);
		}
		for(int j=0; j<columns; j++){
			printf("2D x2=%f", x2[j]);
		}

		//  Fill the 2d array with the sum of squares for x1 and x2
		//  as an example of how to populate the resulting array

		for (int i=0; i<rows; i++){
			for(int j=0; j<columns; j++){
				array2d[i][j] = x1[i]*x1[i] + x2[j]*x2[j];
			}
		}

		//  Print out the 2D array just stored

		printf("\nPrintout of full 2d Array f(x1,x2)\n");
		for(int j=0; j<columns; j++){
			for (int i=0; i<rows; i++){
				printf("row %d column %d: %f", i, j, array2d[i][j]);
			}
		}


		//  Call spline2 to set up 2D 2nd derivative array

		si.spline2(x1, x2, array2d);

		//  Print 2nd derivative array
		printf("\nPrintout of 2nd Derivative Array\n");

		for (int i=0; i<rows; i++){
			for (int j=0; j<columns; j++) {
				printf("2nd deriv:  %d %d: %f", i, j, si.getSpline2d2(i, j));
			}
		}
		printf("\n");


		// Now use splint2 to interpolate particular values from 2D array

		double xx1 = 2.2; 
		double xx2=3.2;
		double yinterpolated = si.splint2(xx1, xx2);

		printf("2D Interpolation: x1=%f x2=%f y=%f", xx1, xx2, yinterpolated);

		xx1 = 1.1; xx2=2.1;
		yinterpolated = si.splint2(xx1, xx2);
		printf("2D Interpolation: x1=%f x2=%f y=%f", xx1, xx2, yinterpolated);

		xx1 = 3.6; xx2=4.5;
		yinterpolated = si.splint2(xx1, xx2);
		printf("2D Interpolation: x1=%f x2=%f y=%f", xx1, xx2, yinterpolated);

		xx1 = 2.4; xx2=4.9;
		yinterpolated = si.splint2(xx1, xx2);
		printf("2D Interpolation: x1=%f x2=%f y=%f", xx1, xx2, yinterpolated);

		xx1 = 16.4; xx2=3.1;
		yinterpolated = si.splint2(xx1, xx2);
		printf("2D Interpolation: x1=%f x2=%f y=%f; return -1 due to attempt to interpolate out of table.", xx1, xx2, yinterpolated);

		xx1 = 5; xx2=4;
		yinterpolated = si.splint2(xx1, xx2);
		printf("2D Interpolation: x1=%f x2=%f y=%f", xx1, xx2, yinterpolated);
	}
};

int main(){
	SplineTester st = SplineTester();
	st.test();

	return 0;

}