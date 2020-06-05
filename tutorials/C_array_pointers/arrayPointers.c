
// Reminder on C arrays and pointers, followed by example of using pointers
// to create arrays of GSL vectors.
//
// Compile with (if on Fedora)
//    gcc arrayPointers.c -o arrayPointers  -lgsl -lgslcblas -lm
// and execute then with
//    ./arrayPointers
// Output is to screen.

#include <stdio.h>
#include <gsl/gsl_matrix.h>

int main (void)
{
    // Let's first remind outselves of the relationship between arrays and 
    // pointers in C
    
    // Create an array
    
    double b[5] = {0.0, 1.0, 4.0, 9.0, 16.0};
    
    // Create a pointer of the same type as the array
    
    double *p;
    
    // Set pointer to initial element of array
    
    p = b;    // An array name is a constant pointer to its 1st element
    
    // Print each element of the array using the pointer p
    
    for (int i = 0; i<5; i++){
        printf("\nb[%d]=%5.2f", i, *(p+i));
    }
    printf("\n");
    
    // Alternatively, manipulate using array name as initial address
    
    double sq;
    for (int i = 0; i<5; i++){
        sq = (*(b+i)) * (*(b+i));   // Square of array entry
        printf("\nb[%d]=%5.2f", i, sq);
    }
    printf("\n\n");
    
    
    // Now use the preceding concepts to create an array of GSL 
    // reaction vectors rv[].
    
    int SIZE = 8;          // Number of (vector) elements in array rv[]
    int ISOTOPES = 3;      // Number of components for each vector 
    
    gsl_vector rv[SIZE];   // Array of type gsl_vector to hold GSL vectors
    gsl_vector *rvPt;      // Pointer to rv[] array
    rvPt = rv;             // Set pointer to beginning address of array rv
   
    // Create an array populated with GSL vectors
    
    printf("\nAllocating an array rv[] of GSL vectors\n");
    for(int i=0; i<SIZE; i++){
        
        //Prototype GSL reaction vector
        gsl_vector * v1 = gsl_vector_alloc (ISOTOPES); 
        
        // Set elements of rv[] pointed to by *rvPt equal to GSL vectors
        *(rvPt+i) = *v1;   
    }
    
    // Fill vector component entries created above with some numbers.  
    // To illustrate we fill the components with the sum of the array 
    // index i and vector index j (cast to a double).
    
    printf("\nPopulate vector components of array rv[i]_j with numbers i+j:");
    for (int i = 0; i < SIZE; i++)
    {
        printf("\n\nrv[%d]",i);
        for(int j=0; j<ISOTOPES; j++){
            gsl_vector_set (rvPt+i, j, (double)i + (double)j);
            printf ("\ni=%d j=%d  i+j =%5.2f  rv[%d]_%d =%5.2f", 
                i, j, i, j, (double)i+(double)j, gsl_vector_get (rvPt+i, j));
        }
    }
    
    // Access and print values of vector components populated above for each rv[] entry
    
    printf("\n\nAccess and print each component of array of GSL vectors rv[i]_j:");
    for (int i = 0; i < SIZE; i++)
    {
        printf("\n");
        for(int j=0; j<ISOTOPES; j++){
            
            // Define a pointer that will point to the GSL vector in array entry rv[i].
            
            gsl_vector *vector = rvPt+i;
            
            // Assign the jth component of the vector in rv[i] to a variable
            
            double component = gsl_vector_get (vector, j);
            printf ("\nrv[%d]_%d = %4.2f", i, j, component);
        }
    }
    
    printf("\n\n");
    
    return 0;
}
