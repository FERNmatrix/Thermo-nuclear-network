
// Examples of using GNU Scientific Library (GSL) for handling matrices and vectors in C.

// See example code at
//     https://sites.google.com/a/phys.buruniv.ac.in/numerical/laboratory/example-codes/using-matrices-with-gsl
// and documentation at
//     https://www.gnu.org/software/gsl/doc/html/vectors.html
//     https://www.gnu.org/software/gsl/doc/html/blas.html
// Need the gsl and blas libraries loaded with #include <gsl/gsl_blas.h> and #include <gsl/gsl_matrix.h>.
// Under Fedora linux need to compile with the flags
//    > gcc matrixGSL.c -o matrixGSL  -lgsl -lgslcblas -lm
// to get all required libraries to link.


#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

int main (void)
{
    // Matrix examples
    
    int rows = 2;  // rows
    int cols = 2;  // columns
    int i, j;
    double value;
    
    printf("\nEXAMPLES OF GNU SCIENTIFIC LIBRARY (GSL) AND BASIC LINEAR ALGEBRA SUBROUTINES (BLAS)\n");
    printf("MATRIX AND VECTOR OPERATIONS IN C\n");
    
    // Allocate matrix memory for two matrices
    gsl_matrix * mat1 = gsl_matrix_alloc (rows, cols);
    gsl_matrix * mat2 = gsl_matrix_alloc (rows, cols);

    printf("\n\nEXAMPLES OF GSL MATRIX OPERATIONS\n");
    // Fill matrix mat1 with values
    printf("\nGenerate values of elements for matrix mat1:");
    for (i = 0; i < rows; i++) {
        printf("\nRow %d: ", i);
        for (j = 0; j < cols; j++) {
            value = (double)(i+j)*(i+j);//*(j-2);
            printf("%5.3f ", value);
            gsl_matrix_set (mat1, i, j, value);
        }
    }
    
    // Read values from the matrix mat1
    printf("\n\nRead matrix elements from matrix mat1:\n");
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            value = gsl_matrix_get (mat1, i, j);
            printf ("mat(%d,%d)=%5.2f  ", i, j, value);
        }
        printf("\n");
    }

    
    // ------------------
    // Vector examples
    // ------------------
    
    printf("\n\nEXAMPLES OF GSL VECTOR OPERATIONS\n");
    
    // Define two vectors of length vsize
    int vsize = 10;
    gsl_vector * v1 = gsl_vector_alloc (vsize);
    gsl_vector * v2 = gsl_vector_alloc (vsize);
    
    // Fill the vectors with values for components
    for (i = 0; i < vsize; i++)
    {
        gsl_vector_set (v1, i, (i+1)*3.0);
        gsl_vector_set (v2, i, 1.0/((i+1)*3.0));
    }
    
    // Print values of vector components
    printf("\nComponents of vector v1:");
    for (i = 0; i < vsize; i++)
    {
        printf ("\nv1_%d = %8.5f", i, gsl_vector_get (v1, i));
    }
    printf("\n\nComponents of vector v2:");
    for (i = 0; i < vsize; i++)
    {
        printf ("\nv2_%d = %8.5f", i, gsl_vector_get (v2, i));
    }
    
    // Multiply the vectors v1 and v2; result stored in v1.  Product
    // should be vector with all entries equal to 1 since v1 and v2
    // were constructed so that the components in v2 and the inverse
    // of components in v1.
    
    gsl_vector_mul(v1,v2);
    printf("\n\nComponents of vector product v1.v2:");
    for (i = 0; i < vsize; i++)
    {
        printf ("\nv1.v2_%d = %8.5f", i, gsl_vector_get (v1, i));
    }
    
    // Compare two vectors with integer components
    
    // Allocate vectors rv1 and rv2
    int rsize = 4;
    gsl_vector * rv1 = gsl_vector_alloc (rsize);
    gsl_vector * rv2 = gsl_vector_alloc (rsize);
    gsl_vector * rv2minus = gsl_vector_alloc(rsize);
    
    // Fill components of rv1 and rv2 with integers
    gsl_vector_set (rv1, 0, -1);
    gsl_vector_set (rv1, 1, -1);
    gsl_vector_set (rv1, 2, 1);
    gsl_vector_set (rv1, 3, 0);
    
    gsl_vector_set (rv2, 0, 1);
    gsl_vector_set (rv2, 1, 1);
    gsl_vector_set (rv2, 2, -1);
    gsl_vector_set (rv2, 3, 0);
    
    printf("\n\n\nCompare vectors to see if they are equal up to a sign");
    
    // Write out the components of rv1 and rv2.  Note that gsl_vector_get
    // returns a double, so it needs to be cast as an int to print correctly.
    printf("\n\nComponents of vector rv1:");
    for (i = 0; i < rsize; i++)
    {
        printf ("\nrv1_%d = %d", i, (int)gsl_vector_get (rv1, i));
    }
    printf("\n\nComponents of vector rv2:");
    for (i = 0; i < rsize; i++)
    {
        printf ("\nrv2_%d = %d", i, (int)gsl_vector_get (rv2, i));
    }
    
    // Compare rv1 with rv2, and rv1 with -rv2 to see if vectors
    // are equivalent, or equivalent up to an overall minus sign
    
    // Compare rv1 and rv2
    int k = gsl_vector_equal(rv1,rv2);
    
    // Compare rv1 and -rv2
    gsl_vector_memcpy(rv2minus, rv2);
    gsl_vector_scale(rv2minus,-1);
    int kk = gsl_vector_equal(rv1,rv2minus);
    
    // Print -rv2.  Note that gsl_vector_get returns a double, so it needs 
    // to be cast as an int to print correctly.
    printf("\n\nComponents of -rv2:");
    for (i = 0; i < rsize; i++)
    {
        printf ("\nrv2minus_%d = %d", i, (int)gsl_vector_get (rv2minus, i));
    }
    
    if(k==1){
        printf ("\n\nCOMPARISON: Vectors rv1 and rv2 are equal\n\n");
    } else if(kk==1){
        printf ("\n\nCOMPARISON: Vector rv2 is equal to -rv1\n\n");
    } else {
        printf ("\n\nCOMPARISON: Vectors rv1 and rv2 are not equal\n\n");
    }

    // ---------------------------------------------
    //  Matrix-matrix multiply using level-3 BLAS
    // ---------------------------------------------
    
    // See http://www.if.pwr.edu.pl/~machnik/Useful_links_files/gsl-ref.pdf
    // pp 126-127
    
    printf("\nMATRIX-MATRIX MULTIPLY C = AB USING GSL AND BLAS 3\n");
    
    // 1D Arrays to create matrix entries A, B, C for C = AB
    double a[] = { 2.0, 3.0, 1.0, 2.0 };
    double b[] = { 2.0, 1.0, 2.0, 3.0 };
    double c[] = { 0.0, 0.0, 0.0, 0.0 };
    
    // Create matrix views from arrays: gsl_matrix_view_array(array, rows, columns)
    gsl_matrix_view A = gsl_matrix_view_array(a, 2, 2);
    gsl_matrix_view B = gsl_matrix_view_array(b, 2, 2);
    gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);
    
    printf("\nMatrix A =\n");
    printf ("Row 0:%5.2f %5.2f\n", a[0], a[1]);
    printf ("Row 1:%5.2f %5.2f\n", a[2], a[3]);
    
    printf("\nMatrix B =\n");
    printf ("Row 0:%5.2f %5.2f\n", b[0], b[1]);
    printf ("Row 1:%5.2f %5.2f\n", b[2], b[3]);
    
    // Compute C = A B
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0,
                    &A.matrix, &B.matrix,
                    0.0, &C.matrix);
    
    // Write product matrix elements
    printf("\nProduct matrix C = AB\n");
    printf ("Row 0:%6.2f %6.2f\n", c[0], c[1]);
    printf ("Row 1:%6.2f %6.2f\n", c[2], c[3]);
    
    
    // -----------------------------------------------------
    //  Example:  Matrix-vector multiply using level-2 BLAS
    // -----------------------------------------------------
    
    // See example at https://na-inet.jp/na/gslsample/linear_system.html
    
    int DIM = 2;
    gsl_matrix *M;
    gsl_vector *x, *y;
    
    // Allocate matrix and vector memory
    M = gsl_matrix_alloc(DIM, DIM);
    x = gsl_vector_alloc(DIM);
    y = gsl_vector_alloc(DIM);
    
    // Populate matrix
    gsl_matrix_set(M, 0, 0, 1.0);
    gsl_matrix_set(M, 0, 1, 2.0);
    gsl_matrix_set(M, 1, 0, 3.0);
    gsl_matrix_set(M, 1, 1, 4.0);
    
    // Populate vector
    gsl_vector_set(x, 0, 1.0);
    gsl_vector_set(x, 1, 2.0);
    
    
    // Matrix vector multiply y = M x
    gsl_blas_dgemv(CblasNoTrans, 1.0, M, x, 0.0, y);
    
    printf("\n\nMATRIX-VECTOR MULTIPLY y = Mx USING GSL AND BLAS 2\n");
    
    printf("\n");
    printf("Matrix M (DIM = %d x %d)\n", DIM,DIM);
    for(i = 0; i < DIM; i++)
    {
        printf("Row %d: ", i);
        for(j = 0; j < DIM; j++)
            printf("%5.2f ", gsl_matrix_get(M, i, j));
        printf("\n");
    }
    
    printf("\n");
    printf("Vector x (DIM = %d)\n", DIM);
    for(i = 0; i < DIM; i++)
    {
        printf("Row %d: ", i);
        printf("%5.2f ", gsl_vector_get(x, i));
        printf("\n");
    }
  
    // Print y vector
    printf("\n");
    printf("Vector y = Mx");
    
    for(i=0; i<DIM; i++){
        printf("\nRow %d: %5.2f", i, gsl_vector_get(y,i));
    }
    printf("\n\n");
    
    // Free memory

    gsl_matrix_free (mat1);
    gsl_matrix_free (mat2);
    gsl_vector_free (v1);
    gsl_vector_free (v2);
    gsl_vector_free (rv1);
    gsl_vector_free (rv2);
    gsl_vector_free (rv2minus);
    gsl_matrix_free (M);
    gsl_vector_free (x);
    gsl_vector_free (y);
    
    return 0;
}
