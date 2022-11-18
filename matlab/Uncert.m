function Uncert()


% Manually change the variables below

dataPoints = 1000; %Plot points (lines in data files)
iso = 16;           % Number of isotopes in network

i=0;                % Loop variable = 0
r = 3;              % No. of calculations being used (for the mean X)

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

            %EXTRA ARRAYS
IsoUnc = zeros(i);
Mean = zeros(i);
DevX = zeros(i);
DevXa = zeros(i);
DevXr = zeros(i);

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
    
    %Fill arrays with log value of X 
    for j=1:iso
        Xlog(i,j) = V(7+j);
        Xrlog(i,j) = U(7+j);
        Xalog(i,j) = W(7+j);
    end 
    Tlog(i,j) = U(1);
    
end

%Mass Fraction array X, should be a matrix [dataPoints,Isotopes)
% Get true value of mass fraction: X = 10^(log(x))
for i=1:dataPoints
    for j=1:iso
        X(i,j) = 10^(Xlog(i,j));
        Xa(i,j) = 10^(Xalog(i,j));
        Xr(i,j) = 10^(Xrlog(i,j));
    end
        T(i) = 10^(Tlog(i));    
    
end

% Calculate mean, deviation from the mean and uncertainty
for i=1:dataPoints
    for j=1:iso
        Mean(i,j) = (X(i,j)+ Xa(i,j)+ Xr(i,j))/(r);
        
        DevX(i,j) = (X(i,j) - Mean(i,j))^2; %Sqaure of deviation
        DevXa(i,j) = (Xa(i,j) - Mean(i,j))^2;
        DevXr(i,j) = (Xr(i,j) - Mean(i,j))^2;
        
        IsoUnc(i,j) = ((DevX(i,j) + DevXa(i,j) + DevXr(i,j))/(6))^(1/2);
    end
end
%(IsoUnc);

%returns maximum value of each column (for any matrix, A)
[M] = max(IsoUnc); % Currently set to show max uncertainty in each iso
C = cat(2,Tlog,IsoUnc); %Adds time to uncertainty matrix

% note Tlog has to be same size as IsoUnc, so there are iso-1 columns of
% all 0s before the log Time is entered in the last column of the Tlog
% matrix.

% For a matrix with 16 isotopes, C will have 32 columns. 16-1 (15) columns
% will be filled with 0s, the 16th will have log time and then the rest
% will be the uncertainty of each mass fraction

%The file "Xerror.m" does a similar calculation but is able to sort out the
%time and mass fractions 1 by 1 and does not have the columns of 0s

writematrix(C,'Uncertainty.txt','delimiter',' ');
type Uncertainty.txt;




fclose(fid);
fclose(fid2);
fclose(fid3);
end
