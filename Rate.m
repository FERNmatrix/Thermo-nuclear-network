function Rate()
        n=0;
        NumReactions = 28;
        T = 10^9;
        T9= T / (1*10^9);     %Temperature in billions of degrees Kelvin
        T913 = T9^(1/3);
        T953 = T9^(5/3);
        rate = zeros(n,1)
        

        
        readfile()    
        
             for i=1:NumReactions
                RP0(i) = reactionParams(1)
                RP1 = P1(i)/T9
                RP2 = P2(i)/T913
                RP3 = P3(i)*T913
                RP4 = P4(i)*T9
                RP5 = P5(i)*T953
                RP6 = P6(i)*log(T9)
        
                rate(i,1) = exp(RP0+RP1+RP2+RP3+RP4+RP5+RP6)
        
                % Rk(n) = exp(P0(n)+(P1(n)/T9) + (P2(n)/T913) + (P3(n)*T913) + (P4(n)*T9) + (P5(n)*T953) + (P6(n)*log(T9)))
             end
end