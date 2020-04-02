function AlterAbundances()
% Using equation D4 in appendix D, and the rates calculated from Reaction
% Rate, calculate the change in abundance of the nuclear species

%call the read in file 
readfile()

if n < NumReactions
    
 nr = NumReactingSpecies(n)
 
 switch(nr)
     case 1 % 1 body reactions
         Flux(n) = Rate(n)*Y(reactant1+n);
     case 2 %2 body reactions
         Flux(n) = Rate(n)*Y(reactant1+n)*Y(reactant2+n);
     case 3 %3 body reactio
         Flux(n) = Rate(n)*Y(reactant1+n)*Y(reactant2+n)*Y(reactant3+n);
 end
 
else
    for j=0:Numreactions-1 %CHECK upper bound on this for loop
        nr = NumReactingSpecies(j);
        
        switch(nr)
            case 1 % 1 body reactions
                Flux(n) = Rate(n)*Y(reactant1+n);
            case 2 %2 body reactions
                 Flux(n) = Rate(n)*Y(reactant1+n)*Y(reactant2+n);
             case 3 %3 body reactio
                 Flux(n) = Rate(n)*Y(reactant1+n)*Y(reactant2+n)*Y(reactant3+n);
        end
    end
end

    






%separate D4 into 3 parts depending on # of reacting species

if NumReactingSpecies == 1
   % Ydot = (Sum on j) Eta(i,j)*Lamda(j)*Y(j)
    else if NumReactingSpecies == 2
        %Ydot = (Sum on j) Eta(i,j)*Lamda(j)*Y(j)+ (Sum on j,k)
        %Eta(i,j,k)*p*Na<jk>*Y(j)*Y(k)

    else if NumReactingSpecies == 3
       %Ydot = (Sum on j) Eta(i,j)*Lamda(j)*Y(j) +
      % (Sum on j,k) Eta(i,j,k)*p*Na<jk>*Y(j)*Y(k) +
      % (sum on j,k,l) Eta(i,j,k,l)*(p^2)*Na^2<jkl>*Y(j)*Y(k)*Y(l)
        end
     end
end

% Possible references in Cuda code LINES 150-177 and maybe 179-202
% LOOP  i s,ines 131-283