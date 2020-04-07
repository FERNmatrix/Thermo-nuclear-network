function rate = ratefun ()

NumReactions = 28;

PN = RateReadFile();

        for i = 1:NumReactions
            P0(i) = PN(i,1);
            P1(i) = PN(i,2);
            P2(i) = PN(i,3);
            P3(i) = PN(i,4);
            P4(i) = PN(i,5);
            P5(i) = PN(i,6);
            P6(i) = PN(i,7);
                    
        end
        %disp(P0)
        %disp(P1)
        %disp(P2)
        %disp(P3)
        %disp(P4)
        %disp(P5)
        %disp(P6)
        
        

%Define variables
i=0;
%NumReactions = 28;
T = 10^9;
T9= T / (1*10^9);     %Temperature in billions of degrees Kelvin
T913 = T9^(1/3);
T953 = T9^(5/3);

%Initialize rate arrays to 0
RP0 = zeros(i);
RP1 = zeros(i);
RP2 = zeros(i);
RP3 = zeros(i);
RP4 = zeros(i);
RP5 = zeros(i);
RP6 = zeros(i);
a = zeros(i);
rate = zeros(i,1);

    for i = 1:NumReactions

    RP0(i) = P0(i);
    RP1(i) = P1(i)/T9;
    RP2(i) = P2(i)/T913;
    RP3(i) = P3(i)*T913;
    RP4(i) = P4(i)*T9;
    RP5(i) = P5(i)*T953;
    RP6(i) = P6(i)*log(T9);

    a(i) = RP0(i) + RP1(i) + RP2(i) + RP3(i) + RP4(i) + RP5(i) + RP6(i);
    
    rate(i) = exp(a(i));

    end %for loop
    %disp(rate(28))
end % function

