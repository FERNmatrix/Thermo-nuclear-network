function readfile (rateLibrary_pp)
    
    %initialize variables

    NumReactions = 28;
    n = 0;
    displayInput = 1;   
    
    %Temporary variables for fscanf issues
    v=0;
    temp = zeros(n);
    Znum = zeros(n);
    
    %Preallocate necessary variables to run faster
                RGclass = zeros(n);
				RGmemberIndex = zeros(n);
				reaclibClass = zeros(n);
				NumReactingSpecies = zeros(n);
				NumProducts = zeros(n);
				isEC = zeros(n);
				isReverseR = zeros(n);
				Prefac = zeros(n);
				Q = zeros(n);
                reactantZ = zeros(n);
                reactantN = zeros(n);
                productZ = zeros(n);
                productN = zeros(n);
                ReactantIndex = zeros(n);
                ProductIndex = zeros(n);
                P0 = zeros(n);
                P1 = zeros(n);
                P2 = zeros(n);
                P3 = zeros(n);
                P4 = zeros(n);
                P5 = zeros(n);
                P6 = zeros(n);
                
                
                
   %Open a file for reading  
   fr = fopen('rateLibrary_pp.data','r');
   
   
%First while loop to loop thrugh all the reactions
while (n < NumReactions)
    subindex = 0;
    
    % Lines 1-2 have set number of input values to be read in
            reactionType = fscanf(fr,'%s',1);
            disp(reactionType);
          
            reactionValues = fscanf(fr,'%d %d %d %d %d %d %d %f %f',9);
            disp(reactionValues);
            
            reactionParams = fscanf(fr,'%f %f %f %f %f %f %f',7);
            disp(reactionParams);
           
            %Reacting Proton NUmber
            for i=1:reactionValues(4)
                Znum(i) = fscanf(fr,'%d',1);
            end
            %disp(Znum);
            
            %Reacting Neutron Number
            for i=1:reactionValues(4)
                Nnum(i) = fscanf(fr,'%d',1);
            end
            %disp(Nnum);
            
            %Product Proton Number
            for i=1:reactionValues(5)
                Zproduct(i) = fscanf(fr,'%d',1);
            end
            %disp(Zproduct);
            
            %Product Neutron Number
            for i=1:reactionValues(5)
                Nproduct(i) = fscanf(fr,'%d',1);
            end
            %disp(Nproduct);
            
            % 
            for i=1:reactionValues(4)
            check1(i) = fscanf(fr,'%d',1);
            end
            %disp(check1);
            
            %
            for i=1:reactionValues(5)
            check2(i) = fscanf(fr,'%d',1);
            end
            %disp(check2);

            while (subindex < 8)
                subindex = subindex +1;
                
                switch (subindex)
           
          case 1 
                %Line 1: reaction title and 9 inputs: string, 7 ints, 2 floats
			n = n+1

                  
				RGclass(n) = reactionValues(1);
				RGmemberIndex(n) = reactionValues(2);
				reaclibClass(n) = reactionValues(3);
				NumReactingSpecies(n) = reactionValues(4);
				NumProducts(n) = reactionValues(5);
				isEC(n) = reactionValues(6);
				isReverseR(n) = reactionValues(7);
				Prefac(n) = reactionValues(8);
				Q(n) = reactionValues(9);
	             
                					
			case 2 
                %Line 2: Reaction Parameters: 7 floats

				P0(n) = reactionParams(1);
				P1(n) = reactionParams(2);
				P2(n) = reactionParams(3);
				P3(n) = reactionParams(4);
				P4(n) = reactionParams(5);
				P5(n) = reactionParams(6);
				P6(n) = reactionParams(7);
				
					
			case 3 
                %Line 3: Proton number in reactants: can take up to 4 ints

                
				for mm=1:NumReactingSpecies(n)
					reactantZ(mm,n) = Znum(mm);
                end
						
			case 4 
                %Line 4: Neutron number in reactants: can take up to 4 ints

                
				for mm=1:NumReactingSpecies(n)
				reactantN(mm,n) = Nnum(mm);
                end    
			
            case 5 
                %Line 5: Proton number in products: can take up to 4 ints

                
				for  mm=1:NumProducts(n)
					productZ(mm,n) = Zproduct(mm);
                end
						
			case 6 
                %Line 6: Nuetron num. in products: can take up to 4 ints

               
				for mm=1:NumProducts(n)
					productN(mm,n) = Nproduct(mm);
                end    
					
			case 7
                %Line 7: Reactant Index number: can take up to 4 ints
                
				for mm=1:NumReactingSpecies(n)
					ReactantIndex(mm,n) = check1(mm);
                end    
						
			case 8
                %Line 8: Product Index number: can take up to 4 ints
                
				for mm=1:NumProducts(n)
					ProductIndex(mm,n) = check2(mm);		
                end	
             end %switch
            end %while subindex
            
    
end %while n< NumReactions
fclose(fr);

%Finding the Rates for each reaction given the reaction parameters
T = 10^9;
T9= T / (1*10^9);     %Temperature in billions of degrees Kelvin
T913 = T9^(1/3);
T953 = T9^(5/3);

Rk = zeros(n,1);

    %while n < NumReactions
    for n=1:NumReactions
        a = P0(n);
        b = P1(n)/T9;
        c = P2(n)/T913;
        d = P3(n)*T913;
        e = P4(n)*T9;
        f = P5(n)*T953;
        g = P6(n)*log(T9);
    
        R = exp(a+b+c+d+e+f+g)
       % R(n)

   % Rk(n) = exp(P0(n)+(P1(n)/T9) + (P2(n)/T913) + (P3(n)*T913) + (P4(n)*T9) + (P5(n)*T953) + (P6(n)*log(T9)))
end



end %function
