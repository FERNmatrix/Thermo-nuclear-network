function Analysis() 

% Manually change the variables below
dataPoints = 1000; %No. of plot output points (lines in data file)
iso = 16;           % isotopes in network

i = 0;              % Loop variable set to 0
            % Time arrays
Tlog = zeros(i);
T = zeros(i);
            %Fast arrays
Xlog = zeros(i);
X = zeros(i);
            % Ref arrays
Xrlog = zeros(i);
Xr = zeros(i);
            % Acc arrays
Xalog = zeros(i);
Xa = zeros(i);
            % Differences
ResXfast = zeros(i);
ResXacc = zeros(i);

% Open data file to be read in
fid = fopen('AlphaHydroplot1.data','r'); % FAST
fid2 = fopen('AlphaHydroRef1.data','r'); % REF
fid3 = fopen('AlphaHydroAcc1.data','r'); % ACC

% read in lines from the file
Columns = 7+iso;

while (i < dataPoints)
    i = i+1;
    %Open %ASY-REF FILE
    U = fscanf(fid2,'%f %f %f',Columns);
    U = U';
    
    %Open PE-FAST FILE
    V = fscanf(fid,'%f %f %f',Columns);
    V = V';
    
    %Open PE-ACC FILE
    W = fscanf(fid3,'%f %f %f',Columns);
    W = W';
    
    for j=1:iso
        Xlog(i,j) = V(7+j);
        Xrlog(i,j) = U(7+j);
        Xalog(i,j) = W(7+j);
    end 
    
    Tlog(i) = U(1);
end

%Mass Fraction array X, should be a matrix [dataPoints,Isotopes)
% Get true values for each X, X = 10^(log(X))
for i=1:dataPoints
    for j=1:iso
        X(i,j) = 10^(Xlog(i,j));
        Xa(i,j) = 10^(Xalog(i,j));
        Xr(i,j) = 10^(Xrlog(i,j));
       
        % Absolute error btwn fast-reference and accurate-reference
        % Abslute error is the residual
        ResXfast(i,j) = (abs(X(i,j) - Xr(i,j)));
        ResXacc(i,j) = (abs(Xa(i,j) - Xr(i,j)));
    end
    
        T(i) = 10^(Tlog(i));
end
% use vpa to go to maximum number of decimals (eliminate rounding)
FAST = vpa(max(ResXfast));
ACC = vpa(max(ResXacc));


FAST = FAST';
ACC = ACC';

%returns maximum value in a column and the line number of the max value
%[Mf,If] = max(diffXfast); 
%[Ma,Ia] = max(diffXacc);

%Sum the rows to find maximum absolute error in the network
for i=1:dataPoints
    Sf(i) = sum(ResXfast(i,:));
    Sa(i) = sum(ResXacc(i,:));
end
Sf =Sf';
Sa = Sa';

%Time of the maximum error 
Tmax = max(Sf);
[Mf, If] = max(Sf);
e1 = [Mf];
t1 = [If];

[Ma,Ia] = max(Sa);
e2 = [Ma];
t2 = [Ia];

fprintf('FAST   time: %f = %d\n',Tlog(t1),t1);
for j=1:iso
    fprintf(' diff[%d]: %f\n',j,vpa(ResXfast(t1,j)));
end
%fprintf('Total error at t1= %f\n',e1);

fprintf('\n');

fprintf('ACC    time: %f = %d\n',Tlog(t2),t2);
for j=1:iso
    fprintf(' diff[%d] : %f\n',j,vpa(ResXacc(t2,j)));
end
%fprintf('Total error at t2 = %f\n',e2);

F =  ResXfast;
A = ResXacc;

writematrix(F,'FastRes.txt','delimiter',' ');
type FastRes.txt;

writematrix(A,'AccRes.txt','delimiter',' ');
type AccRes.txt;

fclose(fid);
fclose(fid2);
fclose(fid3);
end %function

