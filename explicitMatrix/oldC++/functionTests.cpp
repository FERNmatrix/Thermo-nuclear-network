
// To test functions, insert this at the end of main in explicitMatrix.cpp


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
      
      //         // Test of Utilities::stringToChar(string)
      //         printf("TEST Utilities::stringToChar(string) to convert string to Char array:\n");
      //         string ss = "Now is the time";
      //         printf("Char array = %s\n\n", Utilities::stringToChar(ss));
      
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
      
      // Test of eulerUpdate(int i, double FplusSum, double FminusSum, double Y, double dt)
      // and asymptoticUpdate(double Fplus, double Fminus, double Y, double dt)
      
      // Forward Euler
      printf("\nTEST of forward Euler updater\n");
      double fplussum = 801.3;
      double fminussum = 800.0;
      double yy = 0.22;
      double dtt = 0.0001;
      double Yupdate = Integrate::eulerUpdate(0,fplussum, fminussum, yy, dtt);
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
      bool btest = Integrate::checkAsy(fminus, yy);
      if(btest){
          printf("\nIsotope is asymptotic since ck>1\n");
      } else {
          printf("\nIsotope is NOT asymptotic since ck<1\n");
      }
      
      // Test of CPU timer by executing a long, pointless loop
      Utilities::testTimerCPU();
      
  }        // End display of test functions if showFunctionTests==1 */
