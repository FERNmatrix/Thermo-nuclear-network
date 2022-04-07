
/* Example of writing a formatted text (ascii) output file for gnuplot */

 #include <iostream>
 #include <math.h>
 
 using namespace std;
 
 #define PI 3.141592654
 
 int main () {
     
     double x, c, s, t, theta;
     
     // Open a file for output
     FILE * pFile;
     pFile = fopen ("gnufile.data","w");
     
     fprintf(pFile, "#  Deg       Cos       Sin       Tan\n");
     
     for (int i=0; i<73; i++){
         theta = 5.0*(double)i;
         x = theta * PI/180.0;  // Degrees to radians
         c = cos(x);
         s = sin(x);
         t = tan(x);
         fprintf(pFile, "%6.2f %7.3e %7.3e %7.3e\n", theta, c, s, t);
     }

     fclose (pFile);
     return 0;
 }
