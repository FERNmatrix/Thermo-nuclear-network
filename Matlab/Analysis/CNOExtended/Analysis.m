function Analysis() 

% Manually change the variables below

dataPoints = 200;    %No. of plot output points (lines in data file)
iso = 16;           % isotopes in network

i = 0;              % Loop variable set to 0
%__________________________________________________________________%
%                           Storage arrays 
%__________________________________________________________________%
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
%_______________________________________________________________________%

%                 READ IN FILES AND STORE/MODIFY ARRAYS
%_____________________________________________________________________%
% Open data fileS to be read in
fid = fopen('Alphafast1.data','r'); % FAST (or INT)
fid2 = fopen('AlphaRef1.data','r'); % REF
fid3 = fopen('AlphaAcc1.data','r'); % ACC

% read in lines from the file
% 7 Is currently the number of columns before mass fraction columns
% For full strcture/format of the file, look at standard ooutput files
Columns = 7+iso;

%______________________________________________________________________%
% While loop reads the data into variable arrays (U,V,W) 
while (i < dataPoints)
    i = i+1;
    %Open/scan ASY-REF FILE
    U = fscanf(fid2,'%f %f %f',Columns);
    U = U';
    
    %Open/scan PE-FAST FILE
    V = fscanf(fid,'%f %f %f',Columns);
    V = V';
    
    %Open/scan PE-ACC FILE
    W = fscanf(fid3,'%f %f %f',Columns);
    W = W';
 % J-loop only writes the mass fractions into the X arrays
 %(r-ref, a-acc). Note values of data are in log(X)   
    for j=1:iso
        Xlog(i,j) = V(7+j);
        Xrlog(i,j) = U(7+j);
        Xalog(i,j) = W(7+j);
    end 
  % Store the time column in T array, (NOTE Log(t))  
    Tlog(i) = U(1);
end

%________________________________________________________________%
%                   Matrix Manipulation
%________________________________________________________________%

% Mass Fraction array X, should now be a matrix [dataPoints x Isotopes]
% Example if using triple alpha w/ 200 plot outputs array size [200 x 3]

% Array is currently in Log value, so to get true difference,
% get rid of logs: X = 10^(log(X)) for each array (X,Xr,Xa)
for i=1:dataPoints
    for j=1:iso
        X(i,j) = 10^(Xlog(i,j));
        Xa(i,j) = 10^(Xalog(i,j));
        Xr(i,j) = 10^(Xrlog(i,j));
       
% Absolute error (Residuals) btwn fast-ref and acc-ref
        ResXfast(i,j) = (abs(X(i,j) - Xr(i,j)));
        ResXacc(i,j) = (abs(Xa(i,j) - Xr(i,j)));
    end
% Remove Log(t) if desired, may make sense to keep t in log form    
        T(i) = 10^(Tlog(i));
end
%__________________________________________________________________%
%                         Looking at Error values

% 2 ways to use: 
% ---- 1: Look at maximum residual of each isotope to get an idea of accuracy
% on an isotope by isotope basis. Can use ResX---- to plot residuals of
% isotopes over time.
% ---- 2: Using sum of residuals, find the time of maximum error in the
% calculation. This gives idea of accuracy across the network as a whole.
% Can use RMS values to describe accuracy with single value
%__________________________________________________________________%

%                                   Part I:
%__________________________________________________________________%
% Find the maximum error in the network for each isotope
% Use of vpa gives maximum accuracy (like 20+ decimals) this is optional

FAST = vpa(max(ResXfast));
ACC = vpa(max(ResXacc));

% Display max value of the absolute error for each isotope in Fast and Acc

FAST = FAST'
ACC = ACC'

% Returns maximum value in a column and the line number of the max value
% The line number will help relate the error and time at which it occurs

%[Mf,If] = max(ResXfast); 
%[Ma,Ia] = max(ResXacc);

%                                   Part II:
%_________________________________________________________________%
% Sum the rows to find the total error in the network at each ouptut
for i=1:dataPoints
    Sf(i) = sum(ResXfast(i,:));
    Sa(i) = sum(ResXacc(i,:));
end
Sf =Sf';
Sa = Sa';

%                   Finding the time at which max error occurs 
Tmax = max(Sf);     % t = maximum sum of residuals
[Mf, If] = max(Sf); % Show value and line nummber of total error
e1 = [Mf];          % Value of error in fast array (e1)
t1 = [If];          % Line number of max error in fast array (t1)

[Ma,Ia] = max(Sa);
e2 = [Ma];
t2 = [Ia];
%__________________________________________________________________%
%                      Root Mean Square of Error
%__________________________________________________________________%

% In j-loop sum the residuals of each isotope at time w/ line # = t1 or t2
for j=1:iso
   RMSUMf(j) = ResXfast(t1,j)^2;
   
   RMSUMa(j) = ResXacc(t2,j)^2;
end
% Sqrt of that sum
RMSf = (sum(RMSUMf))^(1/2);
RMSa = (sum(RMSUMa))^(1/2);
%___________________________________________________________________%
%                          Print out the errors
%___________________________________________________________________%

fprintf('FAST   Log(t): %f = %d/200 RMS Error = %f \n',Tlog(t1),t1, RMSf);
for j=1:iso
    fprintf(' diff[%d]: %f and  X = %f\n',j,ResXfast(t1,j), X(t1,j));
end
%fprintf('Total error at t1= %f\n',e1);

fprintf('\n');

fprintf('ACC   Log(t): %f = %d/200   RMS Error = %f \n',Tlog(t2),t2, RMSa);
for j=1:iso
    fprintf(' diff[%d] : %f and  X = %f\n',j,ResXacc(t2,j), Xa(t2,j));
end
%fprintf('Total error at t2 = %f\n',e2);

F =  ResXfast;
A = ResXacc;

% Wrtite out data to files

%writematrix(F,'FastRes.txt','delimiter',' ');
%type FastRes.txt;

%writematrix(A,'AccRes.txt','delimiter',' ');
%type AccRes.txt;

fclose(fid);
fclose(fid2);
fclose(fid3);
end %function

